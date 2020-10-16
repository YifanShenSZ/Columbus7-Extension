# Configuration state function (CSF)
The CSF in Columbus7 is specified by `mcdrtin` (`cidrtin`) for MCSCF (MRCI). Manual CSF deletion is possible through modifying the `input walk number (0 to end)` line in `mcdrtin` (`cidrtin`), corresponding to `IWALK` (`IWALKR`) in official documentation

`python3 pick.py` prints the CSFs satisfing the user defined selection rule (specified by modifing function `select`)