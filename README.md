# SmoothingSplines
Time Series Estimation General Additive Models

## Overview
This project was for a final project in Nonparametric Econometric Methods during my masters program at the Toulouse School of Economics. Our prompt was to find any kind of data, explain why splines and regularization methods were suited to that type of data, and fit the data to a smoothing spline. The paper examines the foundations of nonparametric regularization, smoothing splines, and cross validation. I present a theoretical analysis of these topics, as well as how they can be applied to time series data more generally. In the end, I fit a multivariate smoothing spline to hourly traffic count data using the generalized additive model framework in R.

## Data Sources
I found hourly traffic count data at the station id level from Kaggle, [here](https://www.kaggle.com/jboysen/us-traffic-2015). The data contained counts, flow direction, type of road, and characteristics of each station. Since the purpose of this project was time series analysis, I limited the scope of the analysis to one station id.

## Tools Used
This paper was built using Latex, and the coding was done in R. MGCV is the library I used in R for applying generalized additive models.
