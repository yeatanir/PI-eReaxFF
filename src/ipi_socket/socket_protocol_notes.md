# Socket Protocol Notes

Future workflow:

1. receive positions and cell from i-PI;
2. pass them to the existing eReaxFF force provider adapter;
3. return energy, forces, virial, and extras;
4. include extras for charges, ACKS2 diagnostics, status codes, and bond orders
   when available.

