

class definition:
                   # /* Macros for the integer array FLAG      */
    C_B = 0x0000   # /* interior obstacle cells                */
    B_N = 0x0001   # /* obstacle cells adjacent to fluid cells */
    B_S = 0x0002   # /* in the respective direction            */
    B_W = 0x0004   
    B_O = 0x0008 
    B_NW= 0x0005    
    B_SW= 0x0006      
    B_NO= 0x0009    
    B_SO= 0x000a

    C_F = 0x0010   # /* fluid cell */

    C_E = 0x1000   # /* empty celle                 */
    C_N = 0x0800   # /* free surface cells          */
    C_S = 0x0400   # /* adjacent to empty cells     */
    C_W = 0x0200   # /* in the respective direction */
    C_O = 0x0100 
    C_W = 0x0300
    C_NS= 0x0c00
    C_SW= 0x0600
    C_NW= 0x0a00
    C_NO= 0x0900
    C_SO= 0x0500
    C_SWO=0x0700 
    C_NSW=0x0e00
    C_NWO=0x0b00
    C_NSO=0x0d00
    C_NSWO=0x0f00


# /* Macros for POISSON, denoting whether there is an obstacle cell */
# /* adjacent to some direction                                     */
"""
    eps_E !(flag[ii+1, jj] < C_F)
    eps_W !(flag[ii-1, jj] < C_F)
    eps_N !(flag[ii, jj+1] < C_F)
    eps_S !(flag[ii, jj-1] < C_F)
"""