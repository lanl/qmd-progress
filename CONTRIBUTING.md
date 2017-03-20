CONTRIBUTING
=========

# Workflow (Pull Requests)

 Nicolas Bock • 18. March 2017

## Setup of repository

	$ git clone git@github.com:lanl/qmd-progress.git
	$ cd bml
	$ git config user.name “Your Name”
	$ git config user.email “name@somewhere.org”

## Create feature branch

- While working on the pull request it lives on a “feature branch”
- You can name that branch anything you want (in this example we will call it “new_feature”)

	$ git checkout -b new_feature master

- … Make your changes, commit at least once …

	$ git push --set-upstream origin new_feature

- The last step creates the branch on GitHub


## On Github

- Go to [QMD- PROGRESS](https://github.com/lanl/qmd-progress)

- Create a pull request by clicking on the **Compare & pull request ** button

- Add a comment and suggest a reviewer

- Click on the **Create pull request ** to make it effective

## Update pull request

- The pull requested can be updated if changes are necessary

- Add new commits on the “new_feature” branch and

	$  git push

- … or --amend existing commit(s) and

	$ git push force


# Some useful Git commands


  Update the remotes
  ------------------
  $ git remote update --prune         # Update remotes
  $ git rebase upstream/master master # Sync local master with upstream master
  $ git push                          # Sync forked master with local master

  Create new feature branch
  -------------------------
  $ git checkout -b new_feature upstream/master
  $ git log --graph --decorate --all

  Commit changes on feature branch
  --------------------------------
  $ git commit --all

  Create feature branch on forked repository
  ------------------------------------------
  $ git push --set-upstream origin new_feature

  Modify commits on feature branch
  --------------------------------
  $ git commit --all --amend
  $ git push --force

  Once pull request is merged (and the feature branch was deleted on GitHub)
  --------------------------------------------------------------------------
  $ git remote update --prune         # Update remotes
  $ git rebase upstream/master master # Sync local master with upstream master
  $ git push                          # Sync forked master with local master
  $ git branch -d new_feature         # Delete local feature branch

  More helpful commands
  ---------------------
  $ git diff master src/system_mod.F90 > system_mod.F90.patch
  $ git cherry-pick 28ec7eed21e862d41cc7d38ade82e8e4218a504e
