module ereaxff_input_config
  use ereaxff_kinds, only: dp
  use ereaxff_particle_groups, only: particle_group_spec, validate_group_spec
  implicit none
  private

  type, public :: simulation_config
    real(dp) :: dt_fs = 0.1_dp
    integer :: steps = 0
    integer :: seed = 1
    character(len=32) :: integrator = "baoab"
    integer :: output_stride = 1
    integer :: restart_stride = 0
    character(len=512) :: output_dir = ""
  end type simulation_config

  type, public :: force_provider_config
    character(len=32) :: provider_type = "ereaxff_oracle"
    character(len=512) :: v1_path = ""
    character(len=512) :: v2_path = ""
    character(len=512) :: work_root = ""
    character(len=128) :: ffield = "ffield"
    character(len=128) :: geo = "geo"
    character(len=128) :: control = "control"
    character(len=128) :: iopt = "iopt"
    character(len=512) :: smoke_binary = ""
    character(len=512) :: repeat_binary = ""
    character(len=512) :: make_geo_script = ""
    character(len=512) :: timeout_runner = ""
  end type force_provider_config

  type, public :: ring_polymer_config
    character(len=32) :: convention = "standard_ipi"
    integer :: default_beads = 1
    character(len=32) :: propagator = "exact"
  end type ring_polymer_config

  type, public :: charge_coupling_config
    character(len=32) :: mode = "bead_local_acks2"
    real(dp) :: charge_tolerance = 1.0e-8_dp
  end type charge_coupling_config

  type, public :: output_config
    logical :: centroid_xyz = .true.
    logical :: beads_xyz = .true.
    logical :: reaxff_compat = .true.
    logical :: diagnostics = .true.
  end type output_config

  type, public :: input_config
    type(simulation_config) :: simulation
    type(force_provider_config) :: force_provider
    type(ring_polymer_config) :: ring_polymer
    type(charge_coupling_config) :: charge_coupling
    type(output_config) :: output
    type(particle_group_spec), allocatable :: groups(:)
  end type input_config

  public :: read_input_config
  public :: validate_input_config

contains

  subroutine read_input_config(path, config)
    character(len=*), intent(in) :: path
    type(input_config), intent(inout) :: config
    character(len=2048) :: raw_line, line, section, inline_text
    type(particle_group_spec) :: current_group
    logical :: in_group
    integer :: unit, ios, group_count

    call reset_input_config(config)
    section = ""
    in_group = .false.
    group_count = 0

    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "read_input_config: unable to open input file"

    do
      read(unit,'(A)',iostat=ios) raw_line
      if (ios /= 0) exit
      line = strip_comment(raw_line)
      if (len_trim(line) == 0) cycle
      line = adjustl(line)

      if (is_section_header(line)) then
        if (in_group) call append_group(config%groups, group_count, current_group)
        in_group = .false.
        section = section_name(line)
        cycle
      end if

      select case (trim(section))
      case ("simulation")
        call parse_simulation_line(config%simulation, line)
      case ("force_provider")
        call parse_force_provider_line(config%force_provider, line)
      case ("ring_polymer")
        call parse_ring_polymer_line(config%ring_polymer, line)
      case ("charge_coupling")
        call parse_charge_line(config%charge_coupling, line)
      case ("output")
        call parse_output_line(config%output, line)
      case ("groups")
        if (line(1:1) == "-") then
          if (in_group) call append_group(config%groups, group_count, current_group)
          call init_group(current_group)
          in_group = .true.
          inline_text = adjustl(line(2:))
          if (len_trim(inline_text) > 0) call parse_group_line(current_group, inline_text)
        else
          if (.not. in_group) error stop "read_input_config: group entry without list item"
          call parse_group_line(current_group, line)
        end if
      case default
        error stop "read_input_config: unknown or missing section"
      end select
    end do

    if (in_group) call append_group(config%groups, group_count, current_group)
    close(unit)

    call validate_input_config(config)
  end subroutine read_input_config

  subroutine validate_input_config(config)
    type(input_config), intent(in) :: config
    integer :: i

    if (len_trim(config%simulation%integrator) == 0) error stop "validate_input_config: integrator empty"
    if (config%simulation%steps < 0) error stop "validate_input_config: steps must be non-negative"
    if (config%simulation%output_stride < 0) error stop "validate_input_config: output_stride must be non-negative"
    if (config%ring_polymer%default_beads < 1) error stop "validate_input_config: default_beads must be positive"
    if (len_trim(config%force_provider%provider_type) == 0) error stop "validate_input_config: force provider empty"
    if (size(config%groups) < 1) error stop "validate_input_config: no groups provided"

    do i = 1, size(config%groups)
      call validate_group_spec(config%groups(i))
    end do
  end subroutine validate_input_config

  subroutine reset_input_config(config)
    type(input_config), intent(inout) :: config

    config%simulation = simulation_config()
    config%force_provider = force_provider_config()
    config%ring_polymer = ring_polymer_config()
    config%charge_coupling = charge_coupling_config()
    config%output = output_config()
    if (allocated(config%groups)) deallocate(config%groups)
  end subroutine reset_input_config

  subroutine init_group(group)
    type(particle_group_spec), intent(out) :: group

    group%name = ""
    group%selector = ""
    group%dynamics = "classical"
    group%mobile = .true.
    group%beads = 1
    group%thermostat = "none"
    group%temperature_K = 0.0_dp
  end subroutine init_group

  subroutine append_group(groups, count, group)
    type(particle_group_spec), allocatable, intent(inout) :: groups(:)
    integer, intent(inout) :: count
    type(particle_group_spec), intent(in) :: group
    type(particle_group_spec), allocatable :: tmp(:)

    count = count + 1
    allocate(tmp(count))
    if (count > 1) tmp(:count - 1) = groups
    tmp(count) = group
    call move_alloc(tmp, groups)
  end subroutine append_group

  subroutine parse_simulation_line(sim, line)
    type(simulation_config), intent(inout) :: sim
    character(len=*), intent(in) :: line
    character(len=1024) :: key, value

    call split_key_value(line, key, value)
    select case (trim(key))
    case ("dt_fs")
      read(value, *) sim%dt_fs
    case ("steps")
      read(value, *) sim%steps
    case ("seed")
      read(value, *) sim%seed
    case ("integrator")
      sim%integrator = strip_quotes(value)
    case ("output_stride")
      read(value, *) sim%output_stride
    case ("restart_stride")
      read(value, *) sim%restart_stride
    case ("output_dir")
      sim%output_dir = strip_quotes(value)
    case default
      error stop "parse_simulation_line: unsupported key"
    end select
  end subroutine parse_simulation_line

  subroutine parse_force_provider_line(force, line)
    type(force_provider_config), intent(inout) :: force
    character(len=*), intent(in) :: line
    character(len=1024) :: key, value

    call split_key_value(line, key, value)
    select case (trim(key))
    case ("type")
      force%provider_type = strip_quotes(value)
    case ("v1_path")
      force%v1_path = strip_quotes(value)
    case ("v2_path")
      force%v2_path = strip_quotes(value)
    case ("work_root")
      force%work_root = strip_quotes(value)
    case ("ffield")
      force%ffield = strip_quotes(value)
    case ("geo")
      force%geo = strip_quotes(value)
    case ("control")
      force%control = strip_quotes(value)
    case ("iopt")
      force%iopt = strip_quotes(value)
    case ("smoke_binary")
      force%smoke_binary = strip_quotes(value)
    case ("repeat_binary")
      force%repeat_binary = strip_quotes(value)
    case ("make_geo_script")
      force%make_geo_script = strip_quotes(value)
    case ("timeout_runner")
      force%timeout_runner = strip_quotes(value)
    case default
      error stop "parse_force_provider_line: unsupported key"
    end select
  end subroutine parse_force_provider_line

  subroutine parse_ring_polymer_line(rp, line)
    type(ring_polymer_config), intent(inout) :: rp
    character(len=*), intent(in) :: line
    character(len=1024) :: key, value

    call split_key_value(line, key, value)
    select case (trim(key))
    case ("convention")
      rp%convention = strip_quotes(value)
    case ("default_beads")
      read(value, *) rp%default_beads
    case ("propagator")
      rp%propagator = strip_quotes(value)
    case default
      error stop "parse_ring_polymer_line: unsupported key"
    end select
  end subroutine parse_ring_polymer_line

  subroutine parse_charge_line(charge, line)
    type(charge_coupling_config), intent(inout) :: charge
    character(len=*), intent(in) :: line
    character(len=1024) :: key, value

    call split_key_value(line, key, value)
    select case (trim(key))
    case ("mode")
      charge%mode = strip_quotes(value)
    case ("charge_tolerance")
      read(value, *) charge%charge_tolerance
    case default
      error stop "parse_charge_line: unsupported key"
    end select
  end subroutine parse_charge_line

  subroutine parse_output_line(output, line)
    type(output_config), intent(inout) :: output
    character(len=*), intent(in) :: line
    character(len=1024) :: key, value

    call split_key_value(line, key, value)
    select case (trim(key))
    case ("centroid_xyz")
      call parse_logical(value, output%centroid_xyz)
    case ("beads_xyz")
      call parse_logical(value, output%beads_xyz)
    case ("reaxff_compat")
      call parse_logical(value, output%reaxff_compat)
    case ("diagnostics")
      call parse_logical(value, output%diagnostics)
    case default
      error stop "parse_output_line: unsupported key"
    end select
  end subroutine parse_output_line

  subroutine parse_group_line(group, line)
    type(particle_group_spec), intent(inout) :: group
    character(len=*), intent(in) :: line
    character(len=128) :: key, value

    call split_key_value(line, key, value)
    select case (trim(key))
    case ("name")
      group%name = strip_quotes(value)
    case ("selector")
      group%selector = strip_quotes(value)
    case ("dynamics")
      group%dynamics = strip_quotes(value)
    case ("mobile")
      call parse_logical(value, group%mobile)
    case ("beads")
      read(value, *) group%beads
    case ("thermostat")
      group%thermostat = strip_quotes(value)
    case ("temperature_K")
      read(value, *) group%temperature_K
    case ("gamma_fs_inv", "centroid_gamma_fs_inv")
      continue
    case default
      error stop "parse_group_line: unsupported key"
    end select
  end subroutine parse_group_line

  subroutine split_key_value(line, key, value)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key, value
    integer :: colon

    colon = index(line, ":")
    if (colon <= 1) error stop "split_key_value: missing colon"
    key = adjustl(trim(line(:colon - 1)))
    value = adjustl(trim(line(colon + 1:)))
  end subroutine split_key_value

  subroutine parse_logical(text, value)
    character(len=*), intent(in) :: text
    logical, intent(out) :: value
    character(len=:), allocatable :: compact

    compact = to_lower(strip_quotes(text))
    select case (trim(compact))
    case ("true", "t", "yes", "on", "1")
      value = .true.
    case ("false", "f", "no", "off", "0")
      value = .false.
    case default
      error stop "parse_logical: invalid logical value"
    end select
  end subroutine parse_logical

  pure function strip_comment(line) result(out)
    character(len=*), intent(in) :: line
    character(len=:), allocatable :: out
    integer :: hash_pos

    hash_pos = index(line, "#")
    if (hash_pos > 0) then
      out = trim(line(:hash_pos - 1))
    else
      out = trim(line)
    end if
  end function strip_comment

  pure function strip_quotes(text) result(out)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: out
    character(len=:), allocatable :: compact
    integer :: n

    compact = trim(adjustl(text))
    n = len_trim(compact)
    if (n >= 2) then
      if ((compact(1:1) == "'" .and. compact(n:n) == "'") .or. &
          (compact(1:1) == '"' .and. compact(n:n) == '"')) then
        out = compact(2:n - 1)
        return
      end if
    end if
    out = compact
  end function strip_quotes

  pure function to_lower(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, ich

    out = text
    do i = 1, len(text)
      ich = iachar(out(i:i))
      if (ich >= iachar('A') .and. ich <= iachar('Z')) out(i:i) = achar(ich + 32)
    end do
  end function to_lower

  logical function is_section_header(line) result(is_header)
    character(len=*), intent(in) :: line

    is_header = index(line, ":") == len_trim(line)
  end function is_section_header

  pure function section_name(line) result(name)
    character(len=*), intent(in) :: line
    character(len=:), allocatable :: name

    name = to_lower(adjustl(trim(line(:index(line, ":") - 1))))
  end function section_name

end module ereaxff_input_config
