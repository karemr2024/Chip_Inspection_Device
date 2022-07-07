import shutil

source = "C:\Tepegoz\Images"
destination = "C:\Tepegoz\iRiS_Kinetics_Github"

# Copy the content of
# source to destination
dest = shutil.copy(source, destination)

# Print path of newly
# created file
print("Destination path:", dest)

