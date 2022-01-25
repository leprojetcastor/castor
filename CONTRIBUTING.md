# Contributing Code

## How to contribute
The preferred way to contribute to Castor is to fork the [main repository](https://github.com/leprojetcastor/castor) on Github.

1. *Search or open an [issue](https://github.com/UCD4IDS/WaveletsExt.jl/issues) for this project*:   
    * Search for and click into a related issue that you would like to contribute to. If the issue has not been assigned to anyone, leave a comment and assign the issue to yourself by clicking on the 'Assignees' button on the right.

    * If no such issues have been raised, open a new issue by clicking on the 'New issue' on the top right corner. Write a brief description of the issue and wait for the feedbacks from core developers before moving on to the next steps.

1. *Fork the [project repository](https://github.com/UCD4IDS/WaveletsExt.jl)*: Click on the 'Fork' button on the top right corner of the page. This creates a copy of the code under your account on the GitHub server. 

1. *Clone this copy to your local disk*: Open up a terminal on your computer, navigate to your preferred directory, and copy the following.
```
$ git clone git@github.com:<YourGithubUsername>/WaveletsExt.jl.git
$ cd WaveletsExt
```

4. *Instantiate the project*: Open up your Julia REPL, and instantiate the current project in the package manager.
```
$ julia
```

```julia
julia> ]
(v1.7) pkg> activate .
(WaveletsExt) pkg> instantiate
```

5. *Create a branch to hold your changes*:
```
git checkout -b my-feature
```

6. Work on this copy on your computer using your preferred code editor such as VSCode. Make sure you add the [corresponding tests](https://docs.julialang.org/en/v1/stdlib/Test/) in the `test/` directory where appropriate. Use Git to do the version control. When you're done editing, do the following to record your changes in Git:
```
$ git add modified_files
$ git commit
```

7. *Push your changes* to Github with:
```
$ git push -u origin my-feature
```

8. Finally, go to the web page of your fork of the WaveletsExt.jl repository, and *click 'Pull request'* to send your changes to the maintainers for review. This will send an email to the committers.

(If any of the above seems like magic to you, then look up the [Git documentation](https://git-scm.com/doc) on the web.)

### Feature contribution notes
It is recommended to check that your contribution complies with the following rules before submitting a pull request:

- All public methods should have informative docstrings with sample usage presented.

You can also check for common programming errors by testing your code in the package manager mode:
```julia
(WaveletsExt) pkg> test
```

## Filing bugs and feature requests
We use Github issues to track all bugs and feature requests; feel free to [open an issue](https://github.com/UCD4IDS/WaveletsExt.jl/issues) if you have found a bug or wish to see a feature implemented.

It is recommended to check that your issue complies with the following rules before submitting:

* Verify that your issue is not being currently addressed by other [issues](https://github.com/UCD4IDS/WaveletsExt.jl/issues) or [pull requests](https://github.com/UCD4IDS/WaveletsExt.jl/pulls).

* Please ensure all code snippets and error messages are formatted in appropriate code blocks. See [Creating and highlighting code blocks](https://docs.github.com/en/github/writing-on-github/working-with-advanced-formatting/creating-and-highlighting-code-blocks).

* Please include your Julia and WaveletsExt.jl version. This information can be found by running the following code snippet:
```julia
import Pkg
println("Julia $VERSION")       # Julia version
Pkg.status("WaveletsExt")       # Package version
```

## Documentation
You can edit the documentation using any text editor and then generate the HTML output by doing:
```
$ julia --project=docs/ docs/make.jl
```
The resulting HTML files will be placed in `docs/build/` and are viewable in a web browser. See the documentation from [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/man/guide/) for more information.

*Note: Do not commit any files from the `docs/build/` directory. These webpages will be generated automatically via Github Actions when the commit is merged with the master branch of the project repository.*

*Tip: Generally, only the `docs/make.jl` and the files in `docs/src/` need to be updated.*

## Note
This document was gleefully borrowed from [librosa](https://librosa.org/doc/latest/index.html) and [scikit-learn](https://scikit-learn.org/stable/).
