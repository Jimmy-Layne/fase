"""
This function calls our full volume budget for the specified date

This function can be used to run the volume budget for a single day for debugging.

"""
import argparse as arg

import lib.Budget_V3 as bm



# This is a small debug wrapper for a single flight_day
def fase_debug(date):

    # Create budget
    phi_1, zi_1, prf_1, bg_1 = bm.volumeBudget(date)

    # Output data
    output = {'phi':phi_1, 'z':zi_1, 'prf':prf_1, 'bdg':bg_1}
    return output


if __name__ == '__main__':
	# Create Command line arguments
	parser = arg.ArgumentParser(description = "Dates for debugging")
	parser.add_argument("date", type=str, help="Date to be used for debug")
	args = parser.parse_args()

	fase_debug(args.date)


