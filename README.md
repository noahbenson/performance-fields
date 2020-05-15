# Visual Performance Fields

This repository documents the source code and analyses employed in the Visual Performance Fields
project of [Noah C. Benson](https://github.com/noahbenson), [Eline R. Kupers](https://github.com/elinekupers),
Marissa Carrasco, and [Jonathan Winawer](https://github.com/WinawerLab) (2020).

## Contents

The relevant files in this repository fall into two categories: source dode and notebooks.

* The source code is organized as a small python library in the `visual_performance_fields/`
  directory. The directory contains only one file, `__init__.py` which is run on library import.
  When loaded, the library initializes a `neuropythy` dataset containing all data relevant to
  this project. The dataset can be accessed like so:

  ```python
  import neuropythy as ny
  # Force the dataset to initialize.
  import visual_performance_fields

  # Get the dataset object.
  vpf = ny.data['visual_performance_fields']

  # See the help text.
  help(vpl)

  # Get a subject.
  sub = vpl.subject[111312]
  ```

  Note that this repository documents the dataset as it was used in the paper. Beginning with the
  paper's publication, `neuropythy` will contain this dataset by default (i.e., without the
  `visual_performance_fields` package being loaded) and this version of the dataset will be
  maintained with any relevant updates or bug-fixes.
* The analysis notebook `notebooks/performance-fields.ipynb` contains all of the analyses performed
  in the paper and the code used to generate the figures.

## Docker

To run the code in this paper, you may use the `Dockerfile` and `docker-compose.yml` files in
the root of this repository. These files describe a docker-image/virtual-machine in which all
of the analyses included in this repository can be run, assuming your docker installation has
sufficient resources (memory, disk-space, and a working internet connection). Data for the
project will be automatically downloaded from the [Open Science Famework](https://osf.io/5gprz/)
as it is needed in the analyses.

Before running the Docker, you will need to obtain a set of credentials for the Human Connectome
Project's database. Instructions on how to do this are given in the section below. Once you
have obtained these, you can run the docker-image. Make sure that your local port 8888 is free,
then perform the following:

```bash
# In bash:
> git clone https://github.com/noahbenson/visual-performance-fields
...
> cd visual-performance-fields
# Export our HCP credentials; if you have credentials in a file:
> export HCP_CREDENTIALS="`cat ~/.hcp-passwd`"
# If you just have them as a key and secret:
> export HCP_CREDENTIALS="<key>:<secret>"
# Start the jupyter notebook server by bringing up the docker
> docker-compose up
```

This last command should produce a large amount of output as the docker container is built
and started. Once it has started, it should finish by printing a web address that looks
somethin like the following:

```
...
Attaching to visual_performance_fields
visual_performance_fields | Executing the command: jupyter notebook
visual_performance_fields | [I 22:28:43.171 NotebookApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
visual_performance_fields | [I 22:28:44.106 NotebookApp] JupyterLab extension loaded from /opt/conda/lib/python3.7/site-packages/jupyterlab
visual_performance_fields | [I 22:28:44.106 NotebookApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
visual_performance_fields | [I 22:28:44.109 NotebookApp] Serving notebooks from local directory: /home/jovyan
visual_performance_fields | [I 22:28:44.110 NotebookApp] The Jupyter Notebook is running at:
visual_performance_fields | [I 22:28:44.110 NotebookApp] http://(58e2ccd31ba9 or 127.0.0.1):8888/?token=e2f1bd8b37c875799a77198bc240af1b32e1ebc967e04801
visual_performance_fields | [I 22:28:44.110 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
visual_performance_fields | [C 22:28:44.116 NotebookApp]
visual_performance_fields |
visual_performance_fields |     To access the notebook, open this file in a browser:
visual_performance_fields |         file:///home/jovyan/.local/share/jupyter/runtime/nbserver-7-open.html
visual_performance_fields |     Or copy and paste one of these URLs:
visual_performance_fields |         http://(58e2ccd31ba9 or 127.0.0.1):8888/?token=e2f1bd8b37c875799a77198bc240af1b32e1ebc967e04801
```

This final line is telling you how to connect to the notebook server. Basically, copy
everything starting with the "`:8888/`" to the end of the line and paste it into your
browser after "`localhost`", so in this case, you would point your browser to
`localhost:8888/?token=e2f1bd8b37c875799a77198bc240af1b32e1ebc967e04801`. This should
connect you to the notebook server. Click on the `notebooks` directory then on the
`performance-fields.ipynb` file to open the notebook. From there, follow the text and
code in the notebook.


### <a name="credentials"></a> Getting HCP Credentials

Neuropythy uses Amazon's S3 service to obtain structural data from the HCP,
which hosts mosts of its public data there. In order to access these data, you
must obtain a set of S3 credentials from the HCP with access to their S3
buckets. To obtain these credentials, follow these instructions:

1. Point your browser to https://db.humanconnectome.org/ --this should load a
   login page with the title "Access HCP Data Releases"
2. Click the "Register" button in the lower right of the "Create an Account"
   section.
3. Fill out the dialog box that pops up and click "Register". You should get
   a verification email; follow any instructions that it contains.
4. You should now be able to go back to https://db.humanconnectome.org/ and
   log in.
5. Once you have logged in, you should see a page titled "Public Connectome
   Data" with a number of cards below it. The top card should be titled
   "WU-Minn HCP Data - 1200 Subjects" Within this card should be a bunch of
   text describing the dataset and some statistics about it. Near the bottom
   of the card is the word "ACCESS:" followed by a button labeled "Data Use
   Terms Required". Click this button and accept the terms of use that
   appear.
6. The "ACCESS:" tag should now be next to a checkmark and a link labeled
   "Open Access Terms Accepted". Just to the right of this link should be a
   button with the Amazon AWS logo (three yellow cubes) that says something
   about Amazon S3 access. Click this button to bring up the AWS Connection
   Manager dialog; it should have a large button that lets you generate an
   AWS Access Key and Secret. Click this button and follow the instructions
   if any.
7. The access key and secret should look something like this:  
   Key: AKAIG8RT71SWARPYUFUS  
   Secret: TJ/9SJF+AF3J619FA+FAE83+AF3318SXN/K31JB19J4  
   (These are not real credentials).  
   Copy the key and secret and paste them into a file in your home
   directory that you can remember. I recommend using ~/.hcp-passwd, as that
   is the file I will assume you have placed your credentials in during my
   tutorial. When you paste the key and secret into the file, separate them
   by a colon (:) character. For the fake key and secret given above, the
   file would contain the following:  
   AKAIG8RT71SWARPYUFUS:TJ/9SJF+AF3J619FA+FAE83+AF3318SXN/K31JB19J4

For general information on configuring neuropythy, including how to setup the HCP
auto-downloading mechanism, see the [neuropythy configuration wiki
page](https://github.com/noahbenson/neuropythy/wiki/Configuration).


### License 

This README file is part of the visual performance fields girhub repository.

This repository is free software: you can redistribute it and/or Modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.




  
  
  
