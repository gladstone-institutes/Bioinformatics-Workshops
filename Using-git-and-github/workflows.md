# git workflows
Higher-level tips and strategies for more effective usage of git.

## Workflow Models
The appropriate workflow model depends on the size of your team and the type of your project. You have many options, e.g., [Atlassian](https://www.atlassian.com/blog/archives/simple-git-workflow-simple), [Driessen](http://nvie.com/posts/a-successful-git-branching-model/), [GitHub](https://guides.github.com/introduction/flow/) and [GitLab](https://docs.gitlab.com/ee/workflow/gitlab_flow.html), but they tend to be designed for projects like websites that are continuously running and supported by large teams, making them more complex than is needed for a one-off bioinformatics project created by a single individual or a small team. For these types of projects, the GitLab workflow is a good starting point, but you should feel free to customize it to suit your specific needs.

`tldr`
* Use feature branches, e.g., `add-cleaning-step` for changes
* When code in feature branch is solid, merge it into `master` branch
* Always keep your `master` branch in "last known good" state (working, production-ready, deployable)

## What to Commit
Don't commit everything: be deliberate. Never commit files with passwords or API keys. Avoid committing build artifacts, temporary and meta files like `*.swp`, `*.DS_Store`, `*.jar` and `*.pyc` files. Be judicious when deciding whether to commit source data. If you can reference the source and have your scripts download the data if it's not already present, you can reduce the size of your git repo.

To tell git to ignore certain types of files, you can add a `.gitignore_global` file to your home directory and a `.gitignore` file to your project directory. Reference [these templates](https://github.com/github/gitignore) for examples. 

`tldr`
Be sure your home directory has a `.gitignore_global` file like one of these: [Linux](https://github.com/github/gitignore/blob/master/Global/Linux.gitignore), [macOS](https://github.com/github/gitignore/blob/master/Global/macOS.gitignore) or [Windows](https://github.com/github/gitignore/blob/master/Global/Windows.gitignore)

## [Tags](https://git-scm.com/book/en/v2/Git-Basics-Tagging)
Whenever you deploy code or run a publishable job, tag your code using [semantic versioning](http://semver.org/).

## [Rebase](https://git-scm.com/book/en/v2/Git-Branching-Rebasing)

To update a feature branch so the changes it contains are applied on top of the latest from `origin/master`:

```
git fetch
git checkout my-new-feat
git rebase origin/master
```
See [this post](https://blog.algolia.com/master-git-rebase/).
