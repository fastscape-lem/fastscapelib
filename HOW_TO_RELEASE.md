# How to Release Fastscapelib

Some instructions on how to make a new release of Fastscapelib
(work-in-progress).

1. Make sure that your local copy of the ``main`` branch is up-to-date with
   https://github.com/fastscape-lem/fastscapelib.

2. Check and update the release notes in the documentation (include the names of
   contributors). Commit the changes with the message ``prepare for release
   vx.y.z`` (where "x.y.z" should be replaced by the actual version number).

3. Run ``tbump x.y.z --dry-run`` from the project's root directory and check the
   changes that will be made. If everything looks good, run ``tbump z.y.z``. If
   you need to install tbump first, see https://github.com/your-tools/tbump.

4. Push the ``main`` branch with the last commits to
   https://github.com/fastscape-lem/fastscapelib.

5. On the GitHub repository https://github.com/fastscape-lem/fastscapelib,
   go in the "Actions" tab, select the "Python wheels" action and run the
   workflow manually from the ``main`` branch using the button from the GUI.

6. Wait for all running CI workflows (on the ``main`` branch) to succeed.

7. Check the release notes on https://fastscapelib.readthedocs.io/ (once the new
   documentation version has been built and uploaded on ReadTheDocs).

8. If everything looks good, go to
   https://github.com/fastscape-lem/fastscapelib/releases and click on the
   "draft a new release" button. Create a new release ``vx.y.z`` using the tag
   ``x.y.z`` and copy the release notes (summary) from the documentation. This
   should trigger building the Python wheels and upload them to
   https://pypi.org/.

9. After some time, a new pull-request should be opened automatically in the
   https://github.com/conda-forge/fastscapelib-feedstock repository in order to
   create new packages on conda-forge. You can also create a new pull-request
   manually.

10. Run ``tbump u.v.w --only-patch`` with the next version number (development)
    and push the changes.
