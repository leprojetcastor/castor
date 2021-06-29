import sphinx_rtd_theme
#import os

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

def setup(app):
    app.add_css_file("main_stylesheet.css")

extensions = ['breathe', 'sphinx_copybutton']
breathe_projects = { 'castor': '../xml' }
templates_path = ['_templates']
html_static_path = ['_static']
source_suffix = '.rst'
project = 'castor'
copyright = '2020, Ecole Polytechnique, Matthieu Aussal, Marc Bakry, Laurent Series'
author = 'Matthieu Aussal, Marc Bakry, Laurent Series'
html_logo = 'castor_logo.png'

exclude_patterns = []
