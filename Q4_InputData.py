
# simulation settings
POP_SIZE = 2000     # cohort population size
SIM_LENGTH = 50    # length of simulation (years)
ALPHA = 0.05        # significance level for calculating confidence intervals
DISCOUNT = 0     # annual discount rate

ADD_BACKGROUND_MORT = False  # if background mortality should be added
DELTA_T = 1       # years

PSA_ON = False      # if probabilistic sensitivity analysis is on

# transition matrix
TRANS_MATRIX = [
    [0.75, 0.15, 0, 0.1],   # well
    [0, 0.1625, 0.4805, 0.357],  # post-stroke
    [0, 0, 1, 0],   # stroke
    [0, 0, 0, 1]]  # dead


# annual drug costs
Zidovudine_COST = 5000
Lamivudine_COST = 2086.0


