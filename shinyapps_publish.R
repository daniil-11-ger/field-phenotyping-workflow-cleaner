# Deploy the interactive sesame workflow dashboard to shinyapps.io
# 1. Create an account at https://www.shinyapps.io
# 2. In shinyapps.io dashboard, copy your token and secret
# 3. Run the setAccountInfo() command shown in your shinyapps.io account
# 4. Run this script from the project root

install.packages("rsconnect", dependencies = TRUE)
library(rsconnect)

# Example only. Replace with the exact command from your shinyapps.io dashboard:
# rsconnect::setAccountInfo(
#   name = "YOUR_ACCOUNT_NAME",
#   token = "YOUR_TOKEN",
#   secret = "YOUR_SECRET"
# )

rsconnect::deployApp(
  appDir = "shiny_app",
  appName = "sesame-phenotyping-workflow",
  appTitle = "Sesame Phenotyping Workflow"
)
