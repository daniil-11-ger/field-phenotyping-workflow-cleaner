# Run this after completing the workflow and save output for publication reproducibility
sink("outputs/tables/session_info.txt")
print(sessionInfo())
sink()
