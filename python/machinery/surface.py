import numpy as np

/"""
/*------------------------------------------------------------------------*/
/* Initialize particles for free boundary problems to define fluid domain */
/*------------------------------------------------------------------------*/
"""

def particleline *INIT_PARTICLES (int *N,int imax,int jmax,
                                     REAL delx,REAL dely,
                                     int ppc,char *problem,REAL **U,REAL **V)
