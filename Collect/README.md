# Collect

After a batch of Columbus7 single point calculation, to collect data:

    python3 collect.py BatchPath NDirectory NState

Single point information including energy, energy gradient, interstate coupling, transition dipole will be appended

Optionally, you may visualize the potential energy surface by:

    python3 surface.py DataPath

All python scripts support -h or --help argument to show details