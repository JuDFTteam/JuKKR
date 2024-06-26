# How to contribute

We are glad you want to contribute to the Jülich KKR programs.

Here are some important resources:

  * Join our [KKR channel](https://iffchat.fz-juelich.de/signup_user_complete/?id=ascuxqyto3r53pn97gxm1kcu8h) on the iff-chat Mattermost server.
  * Check our [Wiki page](https://iffwiki.fz-juelich.de/kkr/doku.php) and don't be afraid to improve and add things that might be of general interest once you found domething out that is missing there.
  * Take a look at our [online source code documentation](https://kkr.iffgit.fz-juelich.de/jukkr) and the [readme page](https://kkr.iffgit.fz-juelich.de/jukkr/page/index.html) given there.
  * Bugs? Check our [issues page](https://iffgit.fz-juelich.de/kkr/kkrjm/issues) to see if the bug is already known. If it is a new one please [submit a new issue](https://iffgit.fz-juelich.de/kkr/jukkr/issues/new?issue%5Bassignee_id%5D=&issue%5Bmilestone_id%5D=).


## Structure of the repository

The *JuKKR* repository is structured in the following way:
  * `source`: source code, structured in subdirectories
  * `cmake`: files containing the cmake settings
  * `tests`: everything needed to run auto-tests
  * `docs`: source code documentation
  * `utils`: collection of useful utilities, other tools can also be found in the [masci-tools](https://github.com/JuDFTteam/masci-tools) repository.
  * `ElementDataBase`: Database of starting potentials, needed by *voronoi*

Please keep this structure and do your developments in its own branch.


## Submitting changes

Create a new branch starting from the `develop` branch (or the `<codename>-develop` branch for the different sub-codes) to implement your changes. From there it should be merged into the main `develop` branch before it is merged into the master branch with the next release. 

Please follow our coding conventions below and always write a clear log message for your commits.

Also remember to update the *UNRELEASED* section of the changelog (given in the file `CHANGELOG.md`).


## Testing

We use [Gitlab's continuous integration](https://about.gitlab.com/features/gitlab-ci-cd/) features to automatically run tests which are intended to ensure the correctness of our code's features. 

Please make sure that you create a new test case for each new feature you implement or for something that is not covered yet. 

We can always use more test coverage.

The tests are collected int he `tests` subdirectory. The gitlab settings (efinition of the test etc.) can be found in the files of the `tests/gitlab-ci` directory.


## Coding conventions

Our code is (mostly) written in Fortran 90. Please keep to that standard and check your version of the code with multiple compilers. We try to optimize for readability:

  * We indent using two spaces.
  * We ALWAYS put spaces after list items and method parameters (i.e. type `(1, 2, 3)`, not `(1,2,3)`), and around operators (`x + 1`, not `x+1`).
  * We list everything imported from a different module using `only` statements.
  * We clean our code up before publishing it and document it well using the [FORD syntax](https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation).


Thanks,
*The Jülich KKR Team*
