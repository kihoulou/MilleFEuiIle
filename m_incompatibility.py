from m_parameters import *

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKCYAN = '\033[96m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
        
def override_parameters():
    """Performs a check on compatibility of physical inputs in ``m_parameters`` module before running the code.

    :returns: Parameters that have a conflict between each other.
    """
    exit_code = False
    count = 0

    if (rank == 0):
        print("\nPerforming check for incompatible parameters...")

    if (BC_heat_transfer[0][0] == "radiation" and nonlinear_heat_equation == False):
        if (rank == 0):
            count += 1
            print("\tBecause of " + f"{HEADER}BC_heat_transfer[0][0] == 'radiation'{ENDC},"\
                  + " change " + f"{WARNING}'nonlinear_heat_equation' --> True{ENDC}")
        exit_code = True

    if (BC_heat_transfer[0][0] == "radiation" and init_cond_profile == False):
        if (rank == 0):
            count += 1
            print("\tBecause of " + f"{HEADER}BC_heat_transfer[0][0] == 'radiation'{ENDC},"\
                  + " change " + f"{WARNING}'init_cond_profile' --> True{ENDC}")
        exit_code = True
        
    
    if (exit_code == True):
        if (rank == 0):
            if (count == 1):
                print(f"{FAIL}{BOLD}\nThere is " + str(count) + f" incompatible parameter.{ENDC}\n")
            else:
                print(f"{FAIL}{BOLD}\nThere are " + str(count) + f" incompatible parameters.{ENDC}\n")
        exit()
    else:
        if (rank == 0):
            print("No incompatible parameters were found.\n")
