# How to contribute

The preferred way to contribute to Castor is to fork the [main repository](https://github.com/leprojetcastor/castor) on Github.

1. *Search or open an [issue](https://github.com/leprojetcastor/castor/issues) for this project*:   
    * Search for and click into a related issue that you would like to contribute to. If the issue has not been assigned to anyone, leave a comment and assign the issue to yourself by clicking on the 'Assignees' button on the right.

    * If no such issues have been raised, open a new issue by clicking on the 'New issue' on the top right corner. Write a brief description of the issue and wait for the feedbacks from core developers before moving on to the next steps.

2. *Fork the [project repository](https://github.com/leprojetcastor/castor)*: Click on the 'Fork' button on the top right corner of the page. This creates a copy of the code under your account on the GitHub server. 

3. *Clone this copy to your local disk*: Open up a terminal on your computer, navigate to your preferred directory, and copy the following.
```
$ git clone git@github.com:<YourGithubUsername>/castor.git
$ cd castor
```

4. *Instantiate the project*: Install Castor and dependencies and compile the demos following the [user guide](https://leprojetcastor.gitlab.labos.polytechnique.fr/castor/installation.html).

5. *Create a branch to hold your changes*:
```
git checkout -b my-feature
```

6. Work on this copy on your computer using your preferred code editor such as Xcode, Visual Studio, VSCode, etc. Make sure you add the [corresponding tests] in the appropriate `demo_xxx.cpp` file. Use Git to do the version control. When you're done editing, do the following to record your changes in Git:
```
$ git add modified_files
$ git commit
```

7. *Push your changes* to Github with:
```
$ git push -u origin my-feature
```

8. Finally, go to the web page of your fork of the castor repository, and *click 'Pull request'* to send your changes to the maintainers for review. This will send an email to the committers.

If any of the above seems like magic to you, then look up the [Git documentation](https://git-scm.com/doc) on the web. It is recommended to check that your contribution complies with the following rule before submitting a pull request: all public methods should have informative docstrings with sample usage presented, compatible with the documentation process.

