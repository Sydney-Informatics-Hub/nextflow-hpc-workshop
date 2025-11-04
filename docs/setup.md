# Set up your computer

In this workshop, we will be using [Pawsey's Setonix HPC](https://pawsey.org.au/systems/setonix/) and [NCI's Gadi HPC](https://nci.org.au/our-systems/hpc-systems). 

The requirements for this workshop are a personal computer with:

- Visual Studio Code (VSCode)
- A web browser

Below, you will find instructions on how to set up VSCode and connect to the HPC system to which you've been assigned.

Each participant will be provided with their training account and password prior to the workshop.
Before the workshop, you must have the following:

1. VSCode installed
2. The necessary VSCode extensions installed
3. Be able to connect to your assigned HPC.

!!! info

    If you require assistance with the setup, please write in the discussion board on the [Google document]().

## Installing Visual Studio Code

Visual Studio Code (VSCode) is a versatile code editor that we will use for the
workshop. We will use VSCode to connect to the VM, navigate the directories,
edit, view and download files.

1. Download VSCode by following the [installation instructions](https://code.visualstudio.com/docs/setup/setup-overview) for your local Operating System.
2. Open VSCode to confirm it was installed correctly.

![](img/vscode_0.png)

## Installing the VSCode extensions

Specific VSCode extensions are required to connect to the VM and make working with Nextflow files easier (i.e. syntax highlighting).

1. In the VSCode sidebar on the left, click on the extensions button (four blocks)

![](img/vscode_extensions.png)

2. In the Extensions Marketplace search bar, search for `remote ssh`. Select **"Remote - SSH"**

![](img/vscode_ssh_extension.png)

3. Click on the blue `Install` button.

![](img/vscode_ssh_install.png)

4. Once installed, you should see a blue bar in the bottom left corner of the screen. This means that the SSH extension was successfully installed.

![](img/vscode_ssh_installed.png)

5. Close the Extensions tab and sidebar

### Connecting to the HPCs

Ensure you have your training details of your assigned system. 

### Connect to Gadi

TODO

### Connect to Setonix

TODO

### Installing the Nextflow extension

Once you have connected to your assigned HPC, you should also install the Nextflow extension, which provides syntax highlighting and can help identify any potential errors in your code. **Note** that this needs to be done **after** you have connected, as you are installing the extension on the **remote** computer, not your local computer or laptop.

1. Ensure you have connected to your assigned HPC as above. You should see the blue bar in the bottom left corner of the window with the text `SSH: <HPC hostname>`.

![](img/vscode_connected.png)

2. Once again, click on the extensions button in the left sidebar (the icon with four blocks)

![](img/vscode_remote_extensions.png)

3. In the Extensions Marketplace search bar, search for `nextflow` and install the **"nextflow"** estension.

![](img/vscode_nextflow_extension.png)

4. You should see a new icon in the left sidebar that looks like a curvy X shape. This means that the Nextflow extension has been installed correctly. You can now close the extensions tab.

![](img/vscode_nextflow_installed.png)

!!! success

    You have now configured VSCode for the workshop!