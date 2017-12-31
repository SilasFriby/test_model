

#########################################################################
############                 MODEL ASSUMPTIONS               ############
#########################################################################


nStates        <- 8 # number of states is fixed!!
stateCombi     <- c('12', '13', '14', '17', '21', '23', '45', '46', '48', '54', '56') # used for P_ij only! Connections in Markov chain
ageMin         <- 20 # minimum age
ageMax         <- 120  # maximum age
dt  <- 1 / 12 # step length when computing probabilities
dtString <- 'month'
projectionYears <- 60

nStatesTech <- 3
stateCombiTech <- c('12', '13', '21', '23')

nStatesFP <- 3
states_FP <- c(4, 5, 6)
