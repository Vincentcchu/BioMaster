#!/bin/csh
# Script to remove config.yaml from git history

cd /cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/biomaster

# Save the current branch name
set current_branch = `git rev-parse --abbrev-ref HEAD`

# Create a backup branch
set timestamp = `date +"%Y%m%d%H%M%S"`
git branch -m backup-$current_branch-$timestamp

# Checkout a new branch
git checkout --orphan temp_branch

# Copy the template config file to config.yaml
cp config.yaml.template config.yaml

# Add all files
git add .

# Commit
git commit -m "Initial commit with template config.yaml"

# Delete the old branch
git branch -D $current_branch

# Rename the temp branch to the original branch name
git branch -m $current_branch

# Force push to remote
echo "Now run: git push -f origin $current_branch"