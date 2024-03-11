# Piecewise continuous sampling: code and data for American Naturalist submission 

This archive contains all of the code and data used for the manuscript titled, "Piecewise continuous sampling: a method for minimizing bias and sampling effort for estimated metrics of animal behavior". Later, this information will be submitted to Dryad after the initial submission to American naturalist. The goal of this manuscript is to establish piecewise continuous sampling as a flexible method for capturing animal behavior and to give suggestions as to how ethologists can optimally sample their data. This work is based on the observation of ant behavior, but should have applications for other animals as well. 

## Data

The file rawDataTask.csv contains second-by-second measures of the tasks 9 Pogonomyrmex californicus ants were performing over 11,041 seconds. This data is uses for nearly all analyses used in this manuscript, everything from the sample size estimate to evaluating error across various interval numbers. 

The file rawDataActivity.csv is extemely similar to that of the previous file, except the behaviors of the same 9 ants were categorized with a slightly different method. Instead of partioning their behaviors into 9 tasks, we instead categorized their behavior into 3 activity levels. This data is used to validate the results drawn from the previous dataset. 

Finally, agreementData.csv contains the results of an experiment where we evaluated the effect of increasing the number of intervals on observational data. Here, two experimentalists studied the behavior of a single ant for 660 seconds, categorizing their behavior in the same manner as rawDataTask.csv. We then compared the results of these two experiments, and calcuated the degree to which the two agreed (number of seconds tasks were identical for both observers / 660). This 660 second sample would be randomly drawn from a larger 3 hour video, and could be divided into I = 1, 2, 4, 8, 16, and 32 non-overlaping intervals. We did this for 4 colonies. 

## Code

All simulations and analyses were performed in R. For brevity, we compiled all data into a single script (piecewiseContinuousSampling.R) which has been partitioned into different segments which correspond to different parts of the paper. Comments have been included in the script to explain logic. 
