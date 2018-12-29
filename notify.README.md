# __notify__

## __[tl;dr]__
`notify` will run any command, wait for it to complete, and send an email to the user once the command is finished. Particularly useful for long-running programs when run in a `screen`

## __[details]__
When running very long programs, it's nice to not have to think to check on them until they're finished. Also, having an email history of commands run can be a useful reference for future work. `notify` provides an accessible mechanism for keeping track of long-running processes and produces a de facto archive of previous commands.

`notify` can store information about the email server in a configuration file - this will be presented as an option to the user automatically. In addition, it can store information about users, to avoid the user having to enter their email address every time the script is run (though this can be avoided in a variety of other ways, e.g. through aliasing). User information may also be specified on a per-run basis (see usage info).

## __[example usage]__
One requirement of `notify` is that the command being run must be wrapped in quotes – while not required for all commands, failing to use quotes risks breaking the function of the script.

Here is an example using `samtools`,

```
$ notify -e user@email.com "samtools view -b tfoetus.stringtie.sam | samtools sort -o tfoetus.stringtie.sorted.bam"
```
which results in the following email send to 'user@email.com':

subject: 
```
[2017.07.05-13.50][trypsin]: 'samtools view -b tfoetus.stringtie.sam | samtools sort -o tfoetus.stringtie.sorted.bam' completed
```
body:
```
Command-line arguments: samtools view -b tfoetus.stringtie.sam | samtools sort -o tfoetus.stringtie.sorted.bam
Total runtime: 2.704 hours
Return value: 0
Location: /mnt/trypsin_helix/glarue/u12/protists/tritrichomonas_foetus
```

## __[background]__