import bottle
import os
import sys

# add your project directory to the sys.path
project_home = u'/home/bug/cog'
if project_home not in sys.path:
    sys.path = [project_home] + sys.path

# make sure the default templates directory is known to Bottle
templates_dir = os.path.join(project_home, 'views/')
if templates_dir not in bottle.TEMPLATE_PATH:
    bottle.TEMPLATE_PATH.insert(0, templates_dir)

# import bottle application
from bottle_app import application