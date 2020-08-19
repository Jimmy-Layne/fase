""" 
This script handles the curses interface to out analysis, currently deprecated
"""

# 3rd party modules
import curses as cs

# Fase project modules
import lib.Budget_V3 as bm
import lib.visualizations as vis
import lib.fase_menu as mnu
import lib.flight_info as fd



'''This is the interactive analysis client for the FASE project'''


# begin client
std_scr = cs.initscr()
cs.echo()

# begin menu loop
while True:
    # Quit flag
    quit = False
    # clear the screen and input to start
    std_scr.erase()
    cs.flushinp()
    # Get Input Date
    in_loc = mnu.start_screen_draw(std_scr)
    std_scr.addstr(in_loc[0],in_loc[1], "Your Choice: ")
    ch = std_scr.getch()

    if ch == ord('a'):
        date=fd.Date[0]
        phi,zi,prf,bdg = bm.volumeBudget(date)
    elif ch == ord('b'):
        date=fd.Date[1]
        phi,zi,prf,bdg = bm.volumeBudget(date)
    elif ch == ord('c'):
        date=fd.Date[2]
        phi,zi,prf,bdg = bm.volumeBudget(date)
    elif ch == ord('q'):
        break

    
    # Begin Analysis module
    while not quit:
        std_scr.erase()
        cs.flushinp()
        in_loc = mnu.home_screen_draw(std_scr)
        std_scr.addstr(in_loc[0],in_loc[1], "Your Choice: ")
        ch = std_scr.getch()

        if ch == ord('a'):
            mnu.budgets_select(std_scr,bdg,date)

        elif ch == ord('b'):
            mnu.aux_select(bdg)

        elif ch == ord('c'):
            mnu.ent_select(std_scr,bdg,date)
        elif ch == ord('d'):
            mnu.model_select(bdg)

        elif ch == ord('q'):
            quit = True











