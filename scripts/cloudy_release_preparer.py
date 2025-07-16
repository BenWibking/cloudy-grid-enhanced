import os
import subprocess
import glob
import requests
import tarfile

import subprocess

"""
Gold To Do Version: Making a Cloudy Release step by step instructions

Required packages: doxygen, pdfkit, pdflatex, wkhtmltopdf

    1. get a fresh copy of master
            >> git pull
            >> git checkout master
            >> git fresh origin

    2. Update the copy right year
            >> find ./ -type f -exec sed -i -e 's/1978-2023/1978-2025/g' {} \;
        validate the changes
            >> grep -rnw . -e '1978' | grep -v 2025 | grep -v Percival | grep -v Draine

    3. Merge needed branches to master and make sure these changes are described in
       NewCXX and the release paper.

    4. Switch to the release branch
            >> git switch release

    5. do a fast forward merge from master onto the release branch with preference to
       incoming changes if there are any merge conflicts.

            >> git merge master --ff-only
        If merge fails with conflicts:
            >> git merge master
            >> git merge -X theirs master
            >> git push origin release
    The merge_master_to_release.sh does the steps till here.

    6. Bring up two terminals. You will need one terminal to run release preparer script,
       and another to check outputs, and follow the instructions provided.

    7. Run the release script from the root of the release branch

            >> python scripts/cloudy_release_preparer.py
        The script will ask if tsuite has been run, enter "n" to run the tsuite. Then
        once the tsuite is has run clean, come back, re-run the script and enter "y".

    8. Create pdfs from tsuite/auto/doc_suite.htm and tsuite/slow/doc_tsuite.htm

    9. Copy hazy1.pdf, hazy2.pdf, hazy3.pdf, and QuickStart.pdf to top of docs directory

    10. Clean tsuite and source:
            >> cd source
            >> make dist clean

            >> cd ../tsuite
            >> ./clean_tsuite.pl

    11. Commit the changes

    12. Copy Doxygen to the data area
       NOTE: These should not go in the tarball

       First make a new subdirectory under doxygen for the new release
            >> ssh <user-name>@nublado.org
            >> cd /var/www/webapps/data_area/doxygen/
            >> mkdir c25.00
       Exit nublado.org
            >> [CTRL d]
       Copy the doxygen tree to the new subdirectory created in nublado
            >> rsync -a doxygen/html/ cmgu228@nublado.org:/var/www/webapps/data_area/doxygen/c25.00/

    13. Copy the release tarball to nublado
       (this script creates one automatically once all directories have been prepped sucessfully)
            >> rsync -avz c25.00.tar.gz <user-name>@nublado.org:/var/www/webapps/data_area/cloudy_releases/c25/

    14. Tag the latest release branch commit
"""

def prep_source():
    os.chdir("./source/")

    current_dir = os.getcwd()
    print("Entered", current_dir)

    # This checks on any header files that are not used.
    command_args = ["./list_headers.pl"]
    print(f"\n Running source/{command_args[0][1:]}")
    subprocess.run(command_args)
    header_summary = glob.glob(f"{current_dir}/headers.txt")[0]
    header_file_list = glob.glob(f"{current_dir}/listfiles.list")[0]
    print(f" Written to \n{header_summary} \n{header_file_list}")

    # This prunes duplicate headers from source files.
    command_args = ["./uninclude-headers.pl"]
    print(f"\n Pruning delicate headers from source files with {command_args[0]}")
    subprocess.run(command_args)

    # This extracts all atomic data references.
    command_args = ["./doc_atomic_data.pl"]
    print("\n Running ", command_args[0])
    subprocess.run(command_args)

    cloudy_executable = glob.glob(f"{current_dir}/cloudy.exe")
    if f"{current_dir}/cloudy.exe" not in cloudy_executable:
        num_cpus = os.cpu_count()
        command_args = ["make", "-j", f"{num_cpus}"]
        print("Making Cloudy executable for later use.")
        subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("\nSource directory ready for release.\n")
    return 0

def prep_doxygen(cloudy_release):
    os.chdir("./doxygen/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # This checks if doxygen is installed, be looking for its version
    try:
        result = subprocess.run(['doxygen', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            result = {"message": f"Doxygen is installed. Version: {result.stdout.strip()}", "x": 0}
        else:
            result = {"message": f"Doxygen is not installed or not found in the system path.\n Error: {result.stderr.strip()}", "x": -1}
    except:
        result = {"message": f"Doxygen is not installed or not found in the system path.", "x": -1}

    print(result["message"])
    
    # If error returned on doxygen version check, attempt to download doxygen
    if result["x"] == -1:
        print("Attempting to download doxygen...")
        try:
            # Install doxygen using Homebrew
            subprocess.check_call(['brew', 'install', 'doxygen'])
            print("Doxygen installed!")
            install_success = True
        except:
            print("Failed to install Doxygen using Homebrew.")
            install_success = False

        if install_success == False:
            print("Attempting to use sudo apt-get instead...")
            try:
                subprocess.check_call(['sudo', 'apt-get', 'update'])
                subprocess.check_call(['sudo', 'apt-get', 'install', '-y', 'doxygen'])
                print("Doxygen installed!")
                install_success = True                
            except:# subprocess.CalledProcessError as e:
                print("Failed to install Doxygen using apt-get.")
                install_success = False
        
        if install_success == False:
            print("Please install doxygen manually before continuing. Moving onto next directory.")
            return 1

    # This creates the Doxygen documentation
    command_args = ["doxygen", "Doxyfile"]
    print("\n Running ", command_args[0], command_args[1])
    subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    doxygen_html = glob.glob(f"{current_dir}/html/index.html")
    doxygen_latx = glob.glob(f"{current_dir}/latex/")
    print(doxygen_html)
    # If documentation .html file exists open for viewing
    if f"{current_dir}/html/index.html" in doxygen_html:
        print("Doxygen successfully configured.\n Openning doxygen html\n")
        subprocess.run(["open", doxygen_html[0]])

    # FOLLOWING BIT NEEDS TO BE REVISED: 
    # I think we need to be in nebula
    # make sure to first test the if the directory is there
#    print("Next we need to copy the doxygen tree to Cloudy's data server.")
#    nublado_username = input("Please enter (https://data.nublado.org/) username: ")
#    #nublado_password = input("password: ")
#    #print(nublado_username, nublado_password)
#    cloudy_data_server = f"{nublado_username}@nublado.org:/var/www/webapps/data_area/doxygen/{cloudy_release}/"
#    command_args = ["rsync", "-avz", "html/", cloudy_data_server]
#    try:
#        print(f"Copying doxygen directory try to @nublado.org:/var/www/webapps/data_area/doxygen/{cloudy_release}")
#        subprocess.run(command_args)
#    except:
#        print("Could not copy doxygen directory to @nublado.org:/var/www/webapps/data_area/doxygen/")

    print("\nDoxygen directory ready for release.\n")
    return 0

def prep_data():
    os.chdir("./data/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    readme_data_file = "README_data.md"
    print(f"\n Please review and update the data/{readme_data_file}.")
    readme_edit_success = input(" Enter \'continue\' when finished review and update, or enter \'error\' to abort release prep: ")
    if readme_edit_success.lower() == "error":
        print("Error encountered, aborting release prep script.")
        return -1

    # This asks user to make sure all compiled data files are up to date.
    command_args = ["./make_data.sh"]
    print(f"\n Running {command_args[0]}, to make sure compiled data files are up to date.")
    subprocess.run(command_args)

    # This asks user to make sure Cloudy citations are up-to-date
    print("\n Please review and update data/citation_cloudy.txt")
    citation_update_success = input(" Enter \'continue\' once Cloudy citations have been updated, otherwise enter \'error\' to abort: ")
    if citation_update_success == "error":
        print("Error: aborting script.")
        return -1
    else:
        # Then test it by creating and running a test with print citation command
        with open("citation_test.in", 'w') as file:
            file.write("test\n")
            file.write("print citation")
        command_args = ["../source/cloudy.exe", "-r", "citation_test"]
        try:
            subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("\n I created and ran a test for the \'print citation\' command.")
            # Ask user to check the output from the test
            print(" Please review data/citation_test.out")
            citation_test_success = input(" Enter \'continue\' if test looks good, otherwise enter \'error\' to abort: ")
            if citation_test_success == "error":
                print("Error: aborting script.")
                return -1
            subprocess.run(["rm", "citation_test*"])
        except:
            print("Error: aborting, something went wrong running Cloudy executable.")
            return -1

    # Ask user to update the citations for the databases used by Cloudy.
    print("\n Please review and update data/citation_data.txt, the file needs to be updated with the latest database versions.")
    citationdata_update_success = input("Enter \'continue\' once Cloudy citations have been updated, otherwise enter \'error\' to abort: ")
    if citationdata_update_success == "error":
        print("Error: aborting script.")
        return -1

    # This makes sure checksums.dat is up to date.
    # If you build Cloudy in one of the sys_xxxx directories you must temporarily
    # copy (or symlink) vh128sum.exe into source.
    vh128sum_executable = glob.glob("../source/vh128sum.exe")
    if vh128sum_executable != []:
        command_args = ["../scripts/generate_checksums.sh"]
        print(f"\n Running {command_args[0]} to update checksums.dat")
        subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # TODO: add test for checksums run success
    else:
        print("Could not find /source/vh128sum.exe. If you build Cloudy in one of the sys_xxxx")
        print("directories, you must temporarily copy (or symlink) vh128sum.exe into source.")
        print("Moving onto next directory.")
        return 1

    print("\nData directory ready for release.\n")
    return 0


def prep_tsuite():
    os.chdir("./tsuite/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # Prepping tsuite/auto
    os.chdir("./auto/")
    print("\n Entered tsuite/auto/")

    # This creates a list of all test cases, including the input commands and a description of its purpose.
    command_args = ["./doc_tsuite.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, to update list of tests in tsuite/auto.")
    subprocess.run(command_args)

    # doc_tsuite.htm is output from above: contains a tab-delimited list of the files.
    # Open it for user to review
    subprocess.run(["open", "./doc_tsuite.htm"])
    print(" Review opened doc_tsuite.htm for tsuite/auto.")
    doc_tsuite_success = input(" Enter \'continue\' if doc_suites look good, otherwise enter \'error\' to abort: ")
    if doc_tsuite_success == "error":
        print("Error encountered, aborting release prep script.")
        return -1

    # This diagnoses results of tsuite run
    command_args = ["./checkall.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, to diagnose results of tsuite run.")
    files_to_check = ["wvlng.txt", "warnings.txt", "close.txt", "debug.txt", "serious.txt", "crashed.txt"]
    not_empty_files = []
    for file in files_to_check:
        if os.path.exists(file):
            if os.path.getsize(file) != 0:
                not_empty_files.append(file)
        else:
            print(f"Warning! {file} not found.")
    if not_empty_files:
        print("Following files are not empty. Please resolve and come back. Moving onto next directory.")
        return 1

    #The first character should be a sharp sign. This script lists all files that do not start with "#".
    # This is an error, and may indicate that the header was not properly produced.
    command_args = ["./CheckPunchSharp.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]} to make sure save files start with a header saying what the column indicates.")
    subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Prepping tsuite/slow with same steps as above
    os.chdir("../slow/")
    print("\n Entered tsuite/slow/")
    command_args = ["./doc_tsuite.pl"]
    print(f"\n Running tsuite/slow/{command_args[0][2:]}, to update list of tests in tsuite/slow.")
    subprocess.run(command_args)

    subprocess.run(["open", "./doc_tsuite.htm"])
    print("Review opened doc_tsuite.htm for tsuite/slow.")
    doc_tsuite_success = input(" Enter \'continue\' if doc_suites look good, otherwise enter \'error\' to abort: ")
    if doc_tsuite_success == "error":
        print("Error encountered, aborting release prep script.")
        return -1

    command_args = ["./checkall.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, to diagnose results.")
    files_to_check = ["wvlng.txt", "warnings.txt", "close.txt", "debug.txt", "serious.txt", "crashed.txt"]
    not_empty_files = []
    for file in files_to_check:
        if os.path.exists(file):
            if os.path.getsize(file) != 0:
                not_empty_files.append(file)
        else:
            print(f"Warning! {file} not found.")
    if not_empty_files:
        print("Following files are not empty. Please resolve and come back. Moving onto next directory.")
        return 1

    command_args = ["./CheckPunchSharp.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]} to make sure save files start with a header saying what the column indicates.")
    #The first character should be a sharp sign. This script lists all files that do not start with "#". This is an error, and may indicate that the header was not properly produced.
    subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    os.chdir("../")

    # Following script has problems, TODO: fix
#    print("\n Creating pdf from new do doc_tsuite.htm files to include in Hazy2.")
#    # Install needed libraries; TODO: following needs to be handled better if subprocess.run does not work in any instance
#    try:
#        import pdfkit
#    except ImportError:
#        subprocess.run(["pip", "install", "pdfkit"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#        import pdfkit
#    
#    try:
#        subprocess.run(["wkhtmltopdf","--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#    except:
#            try:
#                subprocess.run(["brew", "install", "wkhtmltopdf"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#            except:
#                subprocess.run(["sudo", "apt-get", "install", "wkhtmltopdf"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#
#    pdfkit.from_file("./auto/doc_tsuite.htm", "./auto/doc_tsuite.pdf")
#    pdfkit.from_file("./slow/doc_tsuite.htm", "./slow/doc_tsuite.pdf")
#    print(" New doc_tsuite.pdf files created from doc_tsuite.htm files")
#    #TODO: these need to go to Hazy2

    # TODO: Find coverage run, what script does this? 

    os.chdir("./programs/")
    print("\n Entered tsuite/programs/")

    # This will run the set of programs in run_programs.dat with either gcc or icc
    command_args = ["./run_programs.pl"]
    print(f"\n Running tsuite/programs/{command_args[0][2:]}")
    subprocess.run(command_args)
    # TODO: what is needed to be checked here? Add test.

    print("Tsuite directory ready for release.\n")
    return 0


def prep_docs():
    os.chdir("./docs/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # This makes sure LineLabels.txt and SpeciesLabels.txt are up-to-date
    linelable_input_script = "LineLables"
    print(f"\n Running docs/{linelable_input_script}.in")
    subprocess.run(["../source/cloudy.exe", "-r", linelable_input_script])

    os.chdir("./latex/")
    print("\n Entered", "docs/latex/")

    # Ask user to update citation in Hazy
    print("\n Please review and update citation in item \'CloudyReview\' in docs/latex/common/bibliography2.bib.")
    citation_update_success = input(" Enter \'continue\' once citation is updated, otherwise enter \'error\' to abort: ")
    if citation_update_success == "error":
        print("Encoutered error, aborting script.")
        return -1
    
    # Update table atomic data sources.
    print("\n Please review and update table atomic data sources.")
    datasource_update_success = input(" Enter \'continue\' once table atomic data sources are updated, otherwise enter \'error\' to abort: ")
    if datasource_update_success == "error":
        print("Encoutered error, aborting script.")
        return -1

    print("\n Check whether in-press papers in references have appeared. and update Hazy.")
    datasource_update_success = input(" Enter \'continue\' once in-press papers are updated, otherwise enter \'error\' to abort: ")
    if datasource_update_success == "error":
        print("Encoutered error, aborting script.")
        return -1
    
    print("\n Checking for remaining TODOs' in latex")
    subprocess.run(["grep", "-r", "--include=\"*.tex\"", "TODO", "."])
    todos_update = input(" Enter \'continue\' once TODOs are updated, otherwise enter \'error\' to abort: ")
    if todos_update == "error":
            print("Encoutered error, aborting script.")
            return -1

    hazy_pdfs = glob.glob("hazy*.pdf")
    print("\n Hazy pdfs found in docs/latex/: ",hazy_pdfs)
    # This provides an option to skip re-compiliing hazy pdfs if they have already been
    #  compiled in a release script run previously.
    if "hazy1.pdf" in hazy_pdfs and "hazy2.pdf" in hazy_pdfs and "hazy3.pdf" in hazy_pdfs:
        compile_hazy = input("\n Recompile pdfs (y/n)? ")
    else:
        compile_hazy = "y"

    if compile_hazy.lower() == "y":
        print("\n Compiling Hazy latex files.")
        # Check if pdflatex is installed by doing a version check
        try:
            subprocess.run(["pdflatex", "--version"])
        except:
            print("Your system does not have pdflatex. Wait a moment, attempting to install pdflatex...")
            os_sys = input("What OS system are you using (e.g. mac, linux)? ")
            if "mac" in os_sys:
                try:
                    subprocess.run(["brew", "install", "--cask", "mactex"])
                except:
                    print("Unable to install pdflatex, please install and try again.\nMoving onto next directory.")
                    return 1
            if "linux" in os_sys:
                try:
                    subprocess.run(["sudo", "apt", "install", "texlive"])
                except:
                    print("Unable to install pdflatex, please install and try again.\nMoving onto next directory.")
                    return 1

        # This creates pdfs in each hazy subfolder
        command_args = ["./CompileAll.pl"]
        print(f"\n Running docs/latex/{command_args[0][2:]}, to creating Hazy pdf files.")
        subprocess.run(command_args)

    print("Docs directory ready for release.\n")
    return 0


def main():
    print("Before we get started, the full tsuite must be run.")
    tsuite_run = input("Full tsuite has been run? (y/n) Warning: entering \'n\' will start tsuite run. ")
    if tsuite_run.lower() == "n":
        os.chdir("./tsuite/")
        subprocess.run(["./run_parallel.pl"])
        os.chdir("../")
    elif tsuite_run.lower() == "y":
        dir_prep_success = {}
        cloudy_release = input("Enter cloudy release version number (e.g. \'c25.00\'): ")
        dir_prep_success["source"]  = prep_source()
        if dir_prep_success["source"] >= 0: os.chdir("../")
        else: return
        dir_prep_success["doxygen"] = prep_doxygen(cloudy_release)
        if dir_prep_success["doxygen"] >= 0: os.chdir("../")
        else: return
        dir_prep_success["data"]    = prep_data()
        if dir_prep_success["doxygen"] >= 0: os.chdir("../")
        else: return
        dir_prep_success["tsuite"] = prep_tsuite()
        if dir_prep_success["tsuite"] >= 0: os.chdir("../../")
        # TODO: add a prep_scripts() routine
        else: return
        dir_prep_success["docs"]    = prep_docs()
        os.chdir("../")

        print("Summary: \n")
        for dir in dir_prep_success.keys():
            if dir_prep_success[dir] == 0:
                print(f"{dir} prepped successfully.")
            elif dir_prep_success[dir] == -1:
                print(f"{dir} FAILED!")

        if (-1 not in list(dir_prep_success.values())) and (1 not in list(dir_prep_success.values())):
            # Define your parameters
            output_file = f"{cloudy_release}.tar.gz"
            branch_name = "EditReleaseScript" #release
            prefix_dir = f"{cloudy_release}/"

            # Run the git archive command
            print(f"Making {cloudy_release} tarball...")
            subprocess.run([
                "git", "archive",
                "--format=tar.gz",
                f"--output={output_file}",
                f"--prefix={prefix_dir}",
                branch_name
            ], check=True)

            print(f"Tarball created: {output_file}")
        else:
            print("I did not create a release tarball since some directories failed to be prepped.")
    else:
        print("Aborting release prep script!")
        return

if __name__ == "__main__":
    main()