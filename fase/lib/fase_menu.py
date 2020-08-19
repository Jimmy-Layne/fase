import curses as cs
import lib.visualizations as vis
''' This module contains the print routines for the fase menu'''


def start_screen_draw(win):
    
    win.addstr(0,25,"################################################")
    win.addstr(1,25,"#  Physical Mechanisms of Coastal Fog Retreat  #")
    win.addstr(2,25,"#           From Aircraft Observations         #")
    win.addstr(3,25,"#      ------------------------------------    #")
    win.addstr(4,25,"#              Data Analysis Routines          #")
    win.addstr(5,25,"#              Written By Jimmy Layne          #")
    win.addstr(6,25,"################################################")
    win.addstr(9,15,"Please Select The Day to analyze:")
    win.addstr(10,15,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    win.addstr(11,15,"A) - July 22, 2018")
    win.addstr(12,15,"B) - July 26, 2018")
    win.addstr(13,15,"C) - August 4, 2018")
    win.addstr(14,15,"Q) Quit")
    #Location of return cursor
    x_home = 15
    y_home = 17

    return(y_home,x_home)

def home_screen_draw(win):
    
    win.addstr(1,25,"################################################")
    win.addstr(2,50,"   What would you like to do?") 
    win.addstr(3,50,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    win.addstr(4,50,"A) - Budgets")
    win.addstr(5,50,"B) - Auxilliary figures")
    win.addstr(6,50,"C) - Entrainment Analysis")
    win.addstr(7,50,"D) - Model Analysis")
    win.addstr(8,50,"Q) - Quit")
    
    x_home = 50
    y_home = 11

    return(y_home,x_home)


def budgets_select(win,bdg,date):
    quit = False
    
    while not quit:
        win.erase()
        cs.flushinp()
        in_loc = budgets_screen_draw(win)
        win.addstr(in_loc[0],in_loc[1],"Your Choice: ")
        ch = win.getch()
        
        if ch == ord('a'):
            vis.budget(bdg,phi='se',save=False,show=True)
        elif ch == ord('b'):
            vis.budget(bdg,phi='qt',save=False,show=True)
        elif ch == ord('c'):
            quit = True
    

def budgets_screen_draw(win):
    win.addstr(1,25,"################################################")
    win.addstr(2,25,"   Budget Figures")
    win.addstr(3,50,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    win.addstr(4,50,"A) - Moist Static Energy Budget")
    win.addstr(5,50,"B) - Total Water Budget")
    win.addstr(6,50,"C) - go back")

    x_home = 50
    y_home = 9
    return(y_home,x_home)

def ent_select(win,bdg,deb,date):
    quit = False
    c_filter=False
    
    while not quit:
        win.erase()
        cs.flushinp()
        in_loc = ent_analysis_draw(win,c_filter)
        win.addstr(in_loc[0],in_loc[1],"Your Choice: ")
        ch = win.getch()
        
        if ch == ord('a'):
            vis.budget(bdg,cloud_filter=c_filter,phi='se',save=False,show=True)
        elif ch == ord('b'):
            vis.term_comparison(bdg,c_filter)
        elif ch == ord('c'):
            if c_filter:
                vis.entr_scatter(deb['f_efx_w']['se'],deb['f_efx_p']['se'])
            else:
                vis.entr_scatter(deb['efx_w']['se'],deb['efx_p']['se'])
        elif ch == ord('t'):
            if not c_filter:
                c_filter=True
            else:
                c_filter = False
        elif ch == ord('q'):
            quit = True



# dont be hasty....
def ent_analysis_draw(win,cloud_filter):
    win.addstr(1,25,"################################################")
    win.addstr(2,25,"       Entrainment Analysis Menu")
    win.addstr(3,50,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    win.addstr(4,50,"A) - Moist Static Energy Budget")
    win.addstr(5,50,"B) - Term Comparison")
    win.addstr(6,50,"C) - w vs phi scatter plot")
    win.addstr(7,50,"Q) - go back")
    win.addstr(8,50,"-----------------------------------------------")
    if cloud_filter:
        win.addstr(9,50,"in-cloud filter is currently: ON")
    else:
        win.addstr(9,50,"in-cloud filter is currently: OFF")
    win.addstr(10,50,"Press T to toggle cloud filter")
    win.addstr(11,50,"-----------------------------------------------")
    x_home = 50
    y_home = 14

    return(y_home,x_home)

def model_analysis_draw(win):
    win.addstr(1,25,"################################################")
    win.addstr(2,25,"       Model Analysis Menu")
    win.addstr(3,50,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    win.addstr(4,50,"A) - adv only model")
    win.addstr(5,50,"B) - tend only model")
    win.addstr(6,50,"C) - full tendency model")
    win.addstr(7,50,"D) - comparison figure.")
    win.addstr(8,50,"E) - go back")

    x_home = 50
    y_home = 11

    return(y_home,x_home)


def aux_figures_draw(win):
    win.addstr(1,25,"################################################")
    win.addstr(2,25,"       Auxilliary Figures")
    win.addstr(3,50,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    win.addstr(4,50,"A) - Coastline Figure")
    win.addstr(5,50,"B) - Profile")
    win.addstr(6,50,"C) - go back")

    x_home = 50
    y_home = 7

    return(y_home,x_home)

