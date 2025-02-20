# Install pracma if you haven't already
pacman::p_load(pracma, dplyr, readxl, numDeriv)

## datasets of interest 
### struggling to get all species so just templating with !QULO
data <- readxl::read_excel(path = "Datasets/SMBOmeasures/SMBO2_corMetrics_NRMSE.xlsx",
                          sheet = NULL)
# drop extra rows 
data <- data[1:41,]
# Sample data (replace with your actual data)
x <- as.numeric(data$`Buffer size (kilometers)`) # Input (e.g., advertising spend)
y <- as.numeric(data$`Geographic--"Total Buffer" Approach`) # Output (e.g., sales)

# Calculate the first derivative (marginal return) using finite differences
deriv1 <- numeric(length(x))

# Forward difference for the first point
deriv1 <- (y - y) / (x - x)

# Central difference for the middle points
for (i in 2:(length(x) - 1)) {
  deriv1[i] <- (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
}

# Backward difference for the last point
deriv1[length(x)] <- (y[length(x)] - y[length(x) - 1]) / (x[length(x)] - x[length(x) - 1])

# Calculate the second derivative using finite differences on the first derivative
deriv2 <- numeric(length(x))

# Forward difference for the first point
deriv2 <- (deriv1 - deriv1) / (x - x)

# Central difference for the middle points
for (i in 2:(length(x) - 1)) {
  deriv2[i] <- (deriv1[i + 1] - deriv1[i - 1]) / (x[i + 1] - x[i - 1])
}

# Backward difference for the last point
deriv2[length(x)] <- (deriv1[length(x)] - deriv1[length(x) - 1]) / (x[length(x)] - x[length(x) - 1])

# Find where the second derivative changes sign (inflection point)
inflection_points <- x[which(diff(sign(deriv2))!= 0)]

# Identify the point of diminishing returns
# This is often the first inflection point where the second derivative becomes negative
# (i.e., the point where the marginal return starts decreasing)
point_of_diminishing_returns <- inflection_points[which(deriv2[which(x %in% inflection_points)] < 0)]

# Print the result
print(paste("Point of diminishing returns:", point_of_diminishing_returns))

# Create a data frame with the original points and their derivatives
df <- data.frame(x = x, y = y, deriv1 = deriv1, deriv2 = deriv2)

# Print the data frame
print(df)

# Plot the results
plot(df$x, df$y, type = "l", xlab = "Input (x)", ylab = "Output (y)", main = "Point of Diminishing Returns")
points(point_of_diminishing_returns, df$y[df$x == point_of_diminishing_returns], col = "red", pch = 19)
# legend("topleft", legend = c("Output", "Point of Diminishing Returns"), col = c("black", "red"), lty = c(1, NA), pch = c(NA, 19))