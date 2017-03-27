# CONTRIBUTING

We always welcome changes or additions from outside contributors! Please send
us your pull request!

## Workflow (Pull Requests)

## Setup of repository

    $ git clone git@github.com:lanl/qmd-progress.git
    $ cd bml
    $ git config user.name “Your Name”
    $ git config user.email “name@somewhere.org”

## Create feature branch

While working on a pull request it lives on a "feature branch". You can name
that branch anything you want (in this example we will call it "new_feature")

    $ git checkout -b new_feature master

After creating the feature branch make your changes and commit at least once:

    $ git push --set-upstream origin new_feature

(The last step creates the branch on GitHub.)

## On Github

Go to [QMD- PROGRESS](https://github.com/lanl/qmd-progress). Create the pull
request by clicking on the **Compare & pull request** button. You should add a
comment. In addition you can suggest one or several reviewers but this step is
optional. Now click on **Create pull request**.

## Update pull request

The pull request can be updated by either adding more commits to it or by
amending existing commits. For example, after you add new commits on the
“new_feature” branch you can

    $  git push

to update the pull request. If you changed existing commits on the feature
branch by amending or rebasing the commits a force push is necessary to update
the pull request on GitHub:

    $ git push force

## Some useful Git commands

### Update the remotes

    $ git remote update --prune         # Update remotes
    $ git rebase upstream/master master # Sync local master with upstream master
    $ git push                          # Sync forked master with local master

### Create new feature branch

    $ git checkout -b new_feature upstream/master
    $ git log --graph --decorate --all

### Commit changes on feature branch

    $ git commit --all

### Create feature branch on forked repository

    $ git push --set-upstream origin new_feature

### Modify commits on feature branch

    $ git commit --all --amend
    $ git push --force

### Once pull request is merged (and the feature branch was deleted on GitHub)

    $ git remote update --prune         # Update remotes
    $ git rebase upstream/master master # Sync local master with upstream master
    $ git push                          # Sync forked master with local master
    $ git branch -d new_feature         # Delete local feature branch

### More helpful commands

    $ git diff master src/system_mod.F90 > system_mod.F90.patch
    $ git cherry-pick 28ec7eed21e862d41cc7d38ade82e8e4218a504e
