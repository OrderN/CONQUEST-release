## Contributing to CONQUEST

Thank you for your interest in contributing to the development of
CONQUEST.  The [current
roadmap](https://github.com/OrderN/CONQUEST-release/issues) can be
seen on the [Conquest GitHub issues
page](https://github.com/OrderN/CONQUEST-release/issues).  
The features and developments planned are organised through the
[milestones](https://github.com/OrderN/CONQUEST-release/milestones)
If you want to submit a [new
issue](https://github.com/OrderN/CONQUEST-release/issues/new) please
use one of the templates, and provide as much information as possible,
including input files needed to reproduce any unusual or incorrect
behaviour.  We would also be happy for you to
[fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)
the repository to fix a bug, and then submit a pull request.

Please note that, while we will endeavour to fix any bugs, and are
happy to receive suggestions of new features, we are not able to
guarantee any level of support.

### Joining the development team

We welcome new developers.  If you are interested in contributing a
new feature, please suggest it via the [issues
page](https://github.com/OrderN/CONQUEST-release/issues).  If you are
interested in contributing to an existing feature, please contact the
development team.

We follow a
[gitflow](https://nvie.com/posts/a-successful-git-branching-model/)
workflow (also described
[here](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow). 
The key idea is that all developments should be made on a branch
created from the ``develop`` branch.  The overall approach is this:

* We consider the [GitHub
  repository](https://github.com/OrderN/CONQUEST-release) to be the
  central repository (though as Git is distributed, there is no such
  thing as an actual central repository).  We designate this
  repository as ``origin``.
* The ``origin/master`` branch is the main branch which always represents a
  production-ready state.
* The ``origin/develop`` branch is where developments which have been
  delivered for the next release reside
* Feature branches (named ``f-FeatureName``) branch off ``develop``
  and are used to develop new features or functionality for a future
  release (possibly unspecified).  They may reside on a
  developer's local repository only in the early stages, but will need
  to be pushed for collaboration and discussion.  They will merge into
  ``develop``
* Release branches (named ``release-Release``) branch from ``develop``
  and merge into both ``master`` and ``develop``.
  Release branches allow new developments to proceed on ``develop``
  while final changes and bug fixes are made to the release
* Hotfix branches (named ``hotfix-Version``) branch off ``master`` and
  merge into both ``master`` and ``develop``.  They are only used for
  urgent fixes to an existing release.
