# Architecture Notes

The new engine follows an i-PI-like internal separation:

```text
configuration -> system/groups -> PIMD/thermostat/integrator -> force manager -> force provider
```

The eReaxFF adapter is one force provider.  Future i-PI socket support must
reuse that adapter rather than creating an independent force path.

