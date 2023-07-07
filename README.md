# Alaska event reports

Scripts for automatic generation of event files and daily earthquake reports. REMINDER FOR ETHAN: the copy of this repository on the server is the actual live script, remember to test when making changes. 

Executed each day at 01:00 UTC by adding the following line to `/etc/crontab`:
```
\texttt{  0  1  *  *  * efwillia  /home/efwillia/research/earthquakes/live_scripts/monitor_events.py}
```

Note that absolute paths are required throughout to run in the minimal cron shell environment. I have also added relevant parts of my `$PATH` variable:
```
SHELL=/bin/bash
PATH=/sbin:/bin:/usr/sbin:/usr/bin:/home/efwillia/miniconda3/bin:/home/efwillia/miniconda3/condabin:/home/efwillia/.local/bin:/home/efwillia/bin
```

The main script `monitor_events.py` controls the main parameters that might need to be modified, such as the date, active recording directories, and file naming conventions. It calls two functions from `event_tools.py`, which extract the event data and generate the report. The event data are currently being stored in directories named after the origin time of each event, containing a QuakeML file (`event.qml`) and 2-minute HDF5 files for each array's data (e.g. `KKFLS.h5`). The file format and header information should be consistent with the original data. 

The event reports are produced and broadcast by the `generate_report()` function in `event_tools.py` which first calls the shell script `make_report.gmt` and then emails the resulting PDF using the simple Alpine email client (configured with `smtp.uw.edu`). The latter step is likely possible with a Python SMTP library, but Alpine provides a simple way to keep the user's UW login credentials safe through RSA encryption.

### For first time configuration of Alpine, follow these steps:
(1) Install and launch Alpine
```
>>sudo dnf install alpine
>>alpine
```
(2) Navigate to the configuration menu by entering `s` then `c`. Change the name, domain (`uw.edu`), and SMTP server (`smtp.uw.edu`), then exit with `e` and `q`.
(3) Create a password file for Alpine, and launch Alpine using this file.
```
>>touch .pine-passfile
>>alpine -passfile .pine-passfile
```
(4) Send an arbitrary email using `c` to compose and `[CTRL]+x` to send, when prompted enter UW login credentials, let Alpine save these credentials and set a master password. 
(5) Find the password and encrypt it.
```
>>cd .alpine-smime/.pwd
>>mv MasterPassword.key MasterPassword.key.orig
>>openssl rsa -in MasterPassword.key.orig -out MasterPassword.key
```
(6) Now, as long as Alpine is launched pointing to the password file you can send email without being prompted for login credentials. This only works for sending email with SMTP, as most email clients require oauth2 for IMAP which is challenging to configure. Send email from the command line using the following, for example:
```
>>alpine -passfile ~/.pine-passfile -I ^X,y -subject "SUBJ" -attach "ATTCH.pdf" TO@domain.edu < BODY.txt
```


