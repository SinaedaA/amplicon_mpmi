# Command line for beginners
## Introduction: how to become a computer whisperer
Most likely, if you are a regular computer user, you interact with your computer through what is called a Graphical User Interface (GUI), where the mouse lets you click on buttons and menus to get your computer to do something. 

The Command Line Interface (CLI) can be seen as a way to talk to your computer through text, called command-lines (which are always in English, or derived from English). It is literally a translator between you and your computer. 

The first step to talk to your computer is to open a "terminal emulator", the software that will let you talk to your computer. There is usually one or several pre-installed software on most computers, which can include: `gnome-terminal`, `konsole`, `xterm`, ... In this course, we will use `VSCode` (Microsoft), because it includes many functionalities (terminal, text/script editor and file browser). 

Inside the terminal emulator, there are also different "shells" you can use. The first one was developped by Stephen Bourne, and is called "Bourne Shell" (`sh` in short). Later, he upgraded it and made it free and open-source, under the name **B**ourne **A**gain **SH**ell (`bash`). This is also our default shell, but others exist: **zsh**, **pwsh**, **fish**, etc. Which one you use is a matter of preference, but we will stay with bash for the duration of this course. 

## First look at the command line
First, let's look at the default terminal emulator: GNOME Terminal. When you open it this is what you see:

**insert image of a blank terminal** 

The line that you see written down has the following structure: `username@hostname:location$`. The **username** is you (in this case "**INSERT**"), the **hostname** is the system you are logged onto, and the location is where you are in your computer. Finally, "$" delimits the end of the prompt. 

When I say "where you are on your computer", that is because you are now basically a 🐟 swimming inside your folders. The first place you end up in is your home 🏠. This is indicated by the tilde sign `~`, or by `/Users/<name>` (where <name> is replaced with the name of the computer user). In case you are unsure of your home directory, there is actually a pre-defined **global variable** that contains that information in the computer, called `HOME`. To call a variable, you have to precede the variable name (in this case **HOME**) with a `$` sign in the terminal. We will come back to variables a little later.

You can type the following command in your terminal, press <Enter> and see what happens. 

```bash
echo $HOME
## my output: /Users/sinaeda
```

This should show you your `$HOME` directory. 

Now, let's open the VSCode program. To open a terminal, you can use the shortcut "Strg + Shift + `". The terminal should by default open in the lower part of the software window. The upper part is a "Welcome" page, which you can simply quit. It then let's you choose to open a folder, which will be your base folder for your project. I already made the terminal a bit prettier than our previous example, so that now the "command prompt" is shown underneath the directory we are in, and is preceded by a "❯". 

**add picture of VSCode start screen**

## Command syntax logic
Let's dive into what a **command** is, and how it is structured. Usually, it is structure like this:

```bash
command [-argument] [--long-argument] file
```

The `command` is always one word, that will tell the computer what you want to do. `Short arguments` are one-letter arguments that add certain options to the command, and are preceded by "-". You can combine several one-letter arguments into one, and each individual letter will be interpreted separately, because it is preceded by a single "-". On the other hand, `long arguments` are composed of several letters which are interpreted together as a word. Both of these argument types can be referred to as "**flags**". Often, a short flag has a long equivalent as well. Finally, the `file` specifies what you want to use the command on (when the command acts upon a file). 

Notice that in the previous section, I didn't specify a file, because the `echo` command doesn't act on files, but only on "variables". We'll come back to that later. 

Let's try an example, using the `ls` command. As I said earlier, command are derived from the English language, and `ls` is short for "list directory content". If you are ever unsure what a command does, you can check its "manual" using the `man` command, followed by the command you want to know more information about.

```bash
man ls
```

This will show you what `ls` is doing in its most basic form, and all the flags you can use to fine-tune the output into what you really want to see. To exit the manual, you can simply press `q` (for quit). 

```bash
## ls without any flags or filenames will simply list the contents of the directory you are in
# this should be you $HOME right now
ls
```

However, sometimes we want to see more details about the files inside our directory, so we add options/flags. The most commonly used are `-ltrh`: `-l` stands for "long" (list files in long format, see below), `-t` stands for "sort" (sort files by descending time modified, most recently modified first), `-r` stands for "reverse" (reverse the order of the sort, show most recently modified **last**), and `-h` changes the format in which the file size is displayed. 

The **long format** includes 'file mode', 'number of links', 'owner name', 'group name', 'number of bytes in the file', 'abbreviated month', 'day-of-month last modified', 'hour last modifed', 'minute last modified', and the 'pathname'.

```bash
ls -lrth
```

If you want to know the content of another directory (or file), you can add the "path" to that directory (or file) after the command and its flags. 

```bash
ls -lrth ~/Documents/
```

**Important note**: the arguments (or flags) are *case-sensitive*. As you can see on the `man` page for `ls`, it accepts short flags like `-L` and `-l`, which have different effects on the `ls` command. 

## Basic commands
Here are some basic commands that will enable the little 🐟 in you to navigate around in your computer, and move things around, copy files, look and search inside files, etc. 

- `cd`: stands for "change directory". With this you can move around your files. 
- `cat`: stands for "concatenate". You can display the content of a file in the STDOUT (standard output) of your terminal. You can use `cat` on several files simultaneously, by giving it several file names separated by a space. It will then show you the content of these files, in the order that you specified.
- `grep`: stands for "global regular expression print", and is used to find things in a file (a little bit like <Strg + F> in a PDF for example). We'll explore its capabilities a bit later. 
- `wc`: stands for "word count", and enables you to count words inside a file. Using the `-l` flag, you can count lines instead of words. 
- `mkdir`: stands for "make directory". As its name says, you can use to create a new directory at the specified location. 
- `pwd`: stands for "print working directory". This will print the entire path towards where you are now. 
- `clear`: will clear the content of your screen, if it gets too messy. 
- `history`: will give you a history of commands you ran. 
- `cp`: stands for "copy". You can copy file1.txt to file1_copy.txt, simply by running `cp file1.txt file1_copy.txt`. 
- `mv`: stands for "move". This will move a file from one place to another, and also enables you to rename it. For example, if I want to rename my copied file, I can run `mv file1_copy.txt copy_of_file1.txt`. Using `ls` you can see that the first file "disappeared", because it was renamed to "copy_of_file1.txt". In this case I stayed in the same directory, but I can also use it to move the file around, with or without renaming it. For example, `mv copy_of_file1.txt ~/Downloads/` will move the file to the ~/Downloads/ directory, without renaming it. If I wanted to rename it, I can just add the new name after the directory name : `mv copy_of_file1.txt ~/Downloads/new_file1.txt`. 
- `rm`: stands for "remove". This can delete files and directories, given you use the correct flags. As you can imagine, this can be dangerous as it doesn't move the file to the bin, but rather deletes it **permanently** from the computer. In order to make sure we don't delete important stuff, we will change one of the "settings" of our command line, so that, when you use `rm` on a file, it will prompt you automatically to confirm that you really want to delete the file. 

## Bash settings
In your `$HOME` directory, there are files and folders that you can see using `ls`, but there are also hidden files 👻. The name of these files starts with a `.`, which hides them from you, except if you use the `ls -a` command ("a" stands for "all"). Try it out !

There should already be a file called `.bashrc`. If not, we can create it easily in VSCode : `code $HOME/.bashrc`. This will open an Editor in the upper part of the VSCode window, showing you an empty file (if `.bashrc` did not exist yet) or showing you the content of `.bashrc`. 

This file is loaded whenever you open the terminal and start using bash, meaning we can modify default bash settings. For example, we can change the appearance of our terminal (color scheme).

In our case, we are going to create an `alias` for the `rm` command. An **alias** is a way to create new commands, or to change the way bash executes a command. 

Example:
```bash
## Create a new command that will basically call ls:
alias ShowMeMyDirectoryContentPrettyPlease="ls"
## Try executing this ridiculously long command
ShowMeMyDirectoryContentPrettyPlease
## This will give you the same output as ls
# Note: ls still exists as a command and can be used
```

Now of course, we don't really want to replace a short command like `ls` with a long-ass name as in my example. We would like to make it easier or safer (in the case of `rm`) to use certain commands. 

The `rm` command has a flag `-i` which will prompt the user to confirm whether or not they want to delete the file, as a safety measure. We will create an alias that will replace the default execution of `rm` (delete file without asking for confirmation) by `rm -i`. 

You just have to add:

```bash
alias rm="rm -i"
```

To your `.bashrc` file, save it and close it. Now try it out:

```bash
## create temporary file using touch
touch temp_file.txt
## remove the file
rm temp_file.txt
```

_Did it ask for confirmation to delete the file? Why?_ 

As I mentioned, `.bashrc` is executed/loaded when opening a new terminal with the bash interpreter. We did not quit and re-start the terminal, so `.bashrc` wasn't loaded. Instead of quitting and re-opening, we can also load the file with the `source` command:

```bash
source $HOME/.bashrc
```

Now, if you re-try the example above, it should prompt you for confirmation before deleting the file. 

Other common aliases include:

```bash
alias cp='cp -iv'                           # Preferred 'cp' implementation
alias mv='mv -iv'                           # Preferred 'mv' implementation
alias mkdir='mkdir -pv'                     # Preferred 'mkdir' implementation
alias ll='ls -lrth'                         # Preferred 'ls' implementation
```

Try to figure out what they do differently than their basic counterparts !

## Variables and values
According to [Wikipedia](https://en.wikipedia.org/wiki/Variable_(computer_science)#)

## Exercises
### Creating directories

