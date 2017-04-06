# Canadian Bankrupty Rates

Given monthly data from January 1987 to December 2010 on unemployment rate, population, bankruptcy rate, and housing price index in Canada, the goal of this project was to precisely and accurately forecast monthly bankruptcy rates for Canada. The motivation behind accurately forecasting national bankruptcy rates is because it can be of particular interest to national banks, insurance companies, credit-lenders, and politicians.

## Methods

In this project, I experimented using 4 different time series forecasting methods, which include Holt-Winters, SARIMA, SARIMA with exogenous variables (SARIMAX), and VAR. Holt-Winters and SARIMA rely solely on past history of bankruptcy rates, while the other two incorporate other variables into their predictions.

The first step in choosing the optimal model was to split the data from January 1987 to December 2010 into a training set and a testing set. The training set consisted of monthly data from January 1987 to December 2005 and the testing set consisted of monthly data from January 2006 to December 2010.

After training the different models, I calculated and examined the RMSE of the test set. Additionally, each of these models makes several assumptions which should be met for the predictions to be justifiable. The Holt-Winters and VAR methods make assumptions that the residuals have an average of zero, display constant variance, and are uncorrelated with each other over time. The SARIMA and SARIMAX models make these same assumptions, as well as the assumption that the residuals follow a normal distribution. Combining these steps, I searched for a model that minimized the test RMSE and satisfied each of the assumptions made by the model.

In the search for the optimal model, Holt-Winters and VAR methods gave larger test RMSE values and did not satisfy the majority of the model assumptions. The SARIMA method gave moderately small test RMSE values and satisfied model assumptions, but it was the SARIMAX method that gave the lowest test RMSE value and satisfied all the model assumptions.

The final model I chose was SARIMAX(2,1,3)x(2,1,2) with unemployment rate and population as external variables. This resulted in a test RMSE value of 0.002882778.

## File Descriptions

canadian_bankrupty_rates.R displays the code I used in this project in determining the optimal model to forecast monthly bankruptcy rates in Canada.

data.csv contains all the monthly data used in this project.
