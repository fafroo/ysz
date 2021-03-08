using PyPlot

# TESTING DATA
# - (105 ks): 800 - 660 // 1,2,3,4,5            // 0.3 // data_set="HebbWagner_111", data_set="HebbWagner_110",  data_set="HebbWagner_100",
#
# - (?? ks): 800    // 1,2,3,4,5                // 0.3 - 0.5      // "HebbWagner_110_more_biases"  -> worse, but could be somewhat preprocessed
# - (?? ks): 700    // 1,2,3,4,5                // 0.35 - 0.55    // "HebbWagner_110_more_biases"  -> quite bad quality -> need a lot of preprocessing
#
# - (?? ks): 700    // 1,2,3,4,5                // 0.0, 0.1 - 0.55     // "HebbWagner_111_more_biases"
# - (?? ks): 600    // 1,2,3,4,5                // 0.0, 0.1 - 0.55     // "HebbWagner_111_more_biases"
#
# --> (5*3*2 = 24 ks) 700 //  1,2,3,4,5     // 0.35 : 0.05 : 0.55 //  """110, 111, more_biases"""
  
# Structure
# 
# - Saving tool
#   - together with the structure of preprocessed tree
#   - the structure should be custom to user ... using wildcard-style
#   
# - Main Loop (for menu waiting for instructions from the user)
#   -  update the graphics (clear the window and add new graphs - Old + New )
#   -  the new is done with default settins. But the default settings are written down -> easy to be changed (commenting lines, changing values of procedures)
#   -  the changes are written in the terminal? Or do I want to construct new window also with some text field????
#       - NOT YET, but maybe later
#   - enter to apply the new filtering


















