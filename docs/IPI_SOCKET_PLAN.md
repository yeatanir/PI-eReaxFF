# i-PI Socket Plan

Socket support is planned after the standalone engine is validated.  The socket
adapter will receive positions and cell data, call the same eReaxFF oracle
adapter used by standalone runs, and return energy, forces, virial, and extras
such as charges, ACKS2 status, and bond-order diagnostics.

