#!/usr/bin/bash 

: << 'DOCUMENTATION'

=head1 NAME 

pipelinefunk.sh

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 COPYRIGHT AND LICENSE

Copyright Daniel Nilsson, 2010-2011. 

The package is released under the Perl Artistic License.

=head1 DESCRIPTION

pipelinefunk is a lightweight make-like framework for bash. It
provides a few functions to simplify the use of datestamps and
dependencies to determine what analyses need be updated in a
project. This allows for some trivial checkpointing, and for making
minor updates to analyses etc without having to rerun compute heavy
early analyses. The package also has some rudimentary functions for
file tracking, to simplify cleaning and results packaging tasks. See
more documentation for each function below, but beware of the central
caveat in needsUpdate:

=over 4 

=item *

If your run crashes or there are uncaught problems with input data
etc, incomplete or incorrect results files will be produced. These
will then not be rerun upon a restart, unless a dependent data or
script file is updated. So take care to clean any such
incomplete/incorrect files before resuming.

=back 

=head1 APPENDIX

The rest of the documentation details individual functions and quirks.

=cut

DOCUMENTATION

: << 'FUNCTION_DOC'

=head2 needsUpdate(target, prereq [, prereq]*)

Return true (needsupdate=yes) if target does not yet exist, is older
than its prereqs or forceupdate=yes is in effect. set forceupdate=yes
to run all available analyses, even if the file modification times
advise against it.

Note that this will be false if the target exists, but is incomplete,
as in after a previous interruption or crash, uncaught problems with
input data etc.  Take care to remove any affected intermediate or
results files after an error was detected, or clean the whole
directory with directive C<shinyclean> if unsure.

=cut

FUNCTION_DOC

if [ -z "$forceupdate" ]
then
    forceupdate=no
fi

function needsUpdate()
{
    needsupdate="no"
    
    if [ "$forceupdate" = "yes" ] 
    then
	needsupdate="yes"
    fi

    target=$1;
    
    for prereq in ${@:2}
    do
	if [ $target -ot $prereq ]
	then
	    needsupdate="yes"
	fi
   done
    
    [ "$needsupdate" = "yes" ]
}

: << 'NUTB_FUNCTION_DOC'

=head2 needsUpdateTimeBased(target, timetoobsolete)

Return true (needsupdate=yes) if target does not yet exist, is older than timetoobsolete (in seconds) or forceupdate=yes is in effect. set forceupdate=yes to run all available analyses, even if the file modification times advise against it. 

    # sample synopsis
   
    seconds_in_two_days=$(( 60 * 60 * 24 * 2))
    update_pathways=no

    org_kegg_list=${org}.list

    if needsUpdateTimeBased ${org}.list $seconds_in_two_days
    then
	wget -c ftp://ftp.genome.jp/pub/kegg/pathway/organisms/${org}/${org}.list
	update_pathways=yes
	updates=yes
    fi

=cut

NUTB_FUNCTION_DOC

function needsUpdateTimeBased()
{
    local file=$1
    local timetoobsolete=$2
        
    local filestamp
    local nowstamp=`date +%s`

    local needsupdate="no"

    if [ "$forceupdate" = "yes" ] 
    then
	needsupdate=yes
    fi

    if [ ! -w $file ]
    then
	needsupdate=yes
    else
	# stat is used for timestamp retrieval, and works slightly differently on OSX
	os=`uname`
	if [ $os = "Darwin" ]
	then
	# osx
	    filestamp=`stat -f %c $file`
	elif [ $os = "Linux" ] 
	then
	# linux 
	    filestamp=`stat --format %Z $file`
	fi

	if [ $(( $nowstamp - $filestamp )) -gt $timetoobsolete ] 
	then
	    needsupdate=yes
	fi
    fi

    [ "$needsupdate" = "yes" ]
}

: <<'REGISTER_FUNCTION_DOC'

=head2 registerFile(file, category) 

 USAGE: registerFile /absolute/path/to/file category 

Register file with the cleanup system. Basically log it to a file,
named after the category. Category is currently limited to
C<result|temp>. Best practice is to call registerFile before actually
starting to write the file, so that cleanup will find also partial
files.

CAVEAT: registerFile is not inherently multiple-worker safe. There
would be a race condition between inventory of register and addition
of a new item. This is not deemed severe, as the bulk (all?) of the
filenames to be listed will come from work segments that are made
exclusive through the work sharing lock. See workLockOk.

=cut

REGISTER_FUNCTION_DOC

# on load, =)
# save current PWD for use with later regs.
pipelineregisterdir=$PWD

function registerFile()
{
    local savefile=$1
    local category=$2

    register=${pipelineregisterdir}/.pipeline.register.$category

    # create on first use
    if [ ! -e $register ]
    then
	touch "$register"
    fi

    # check that it's not already on the list?
    grep -x "$savefile" ${register} > /dev/null
    if [ $? -eq 0 ]
    then
	# savefile was already on the register file
	:
    else       
	echo "$savefile" >> ${register}
    fi
}

#: <<'FLAG_DOC'
#
#=head2 flagActive(file)
#
#USAGE: flagActive file 
#
#Flags file as being written. Use to provide some protection against inconsistencies when a run is interrupted.
#
#TODO: make a bit more multiple friendly by separating flag list 
#
# actually not very useful without good interruption traps - revisit when such are in place. 
# atomic runs via temp output and moving files after operation successfully completes would be an option, but costs peak disk space
#
#FLAG_DOC
#
#flagfile=${pipelineregisterdir}/.pipeline.$$.active.flag

#function flagActive()
#{
#    local file=$1
#
#    local nowstamp=`date +%s`
#
#    touch ${file}.active.${nowstamp}
#    echo ${file}.active.${nowstamp} > ${flagfile}
#}

#function flagDone()
#{
#    if [ $? != 0 ] ... check return status for a little bit more protection?    
#    rm -f `cat $flagfile`
#    rm -f "$flagfile"
#}

# remove any outstanding "active" files on load.
#if [ -e $flagfile ]
#then
#    rm `cat $flagfile`
#    rm "$flagfile"
#fi

: <<'CLEAN_FUNCTION_DOC'

=head2 cleanCategory(category)

USAGE: cleanCategory category 

Delete files registered with the cleanup system. Will recursively
delete directories if registered. Category is currently limited to
C<result|temp>.

=cut

CLEAN_FUNCTION_DOC

function cleanCategory()
{
    local category=$1

    register=${pipelineregisterdir}/.pipeline.register.${category}
 
    if [ -e "$register" ] 
    then
	for file in `cat $register`
	do
	    if [ -d "$file" ] 
	    then
		rm -rf "$file"
	    else
		rm "$file"
	    fi
	done
	rm "${pipelineregisterdir}/.pipeline.register.${category}"
	echo "Cleaned up ${category} files."
    else
	echo "No register file $register found. Directory was perhaps already clean?"
    fi
}


: <<'WORK_LOCK_OK'

=head2 workLockOk(file)

USAGE: workLockOk file  

Attempt to obtain work lock for file. Return true if lock was
obtained. Return false if lock was not obtained. Use e.g. to enable
multiple workers to co-exist on the same directory. Lock eg raw file
before entering a work segment. Carefully choose the file to lock to
atomically lock all downstream work intended. This avoids both
deadlocks and having another worker start a downstream segment on soon
to be outdated previous results.  

=cut

WORK_LOCK_OK

function workLockOk()
{
    local file=$1
    
    local lockfile=${file}.lock
    local idstring=${HOSTNAME}"-"$$

    if [ -a "$file" ]
    then
	if ( set -o noclobber; echo "$idstring" > "$lockfile") 2> /dev/null;
	then
	    trap 'rm -f "$lockfile"; exit $?' INT TERM EXIT
	    return 0;
	else 
	    return 1;
	fi

    fi
}

: <<'POD_RELEASE_LOCK'

=head2 releaseLock(file)

USAGE: releaseLock file  

Release work lock for file. Only call if an exclusive lock is already in effect.
Will nevertheless perform simple check for ownership (hostname-script exclusive).

=cut

POD_RELEASE_LOCK

function releaseLock()
{
    local file=$1
    
    local lockfile=${file}.lock
    local idstring=${HOSTNAME}"-"$$

    # check that we own the lock (basically this function should not be called without a lock, but better safe...)

    if [ -a "$lockfile" -a "`cat $lockfile`" == "$idstring" ]
    then
 	rm -f "$lockfile"
    fi

    trap - INT TERM EXIT
}

: << 'DOC_DIRECTIVE'

=head1 SYNOPSIS

 [DIRECTIVE='onlyclean|onlyshinyclean|clean|shinyclean']
 . pipelinefunk.sh

If a directive of e.g. C<onlyclean> is set already when sourcing this,
then clean accordingly. A directive can be given in C<$1> or
C<$DIRECTIVE>. Note that for most pipelines, giving the directive on
the command line as C<$1> would interfere with normal usage, and so be
limited to onlyclean and onlyshinyclean.

=over 4

=item C<clean> 

Clean temp files.

=item C<shinyclean>

Clean both temp files and results.

=item C<onlyclean>

Clean temp files, then terminate unconditionally
without returning to the rest of the pipeline.

=item C<onlyshinyclean>

Clean both temp files and results, then terminate unconditionally
without returning to the rest of the pipeline.

=back

=cut

DOC_DIRECTIVE

if [ -z "$DIRECTIVE" ]
then
    DIRECTIVE=$1
fi

if [ "$DIRECTIVE" == "clean" ]
then
    cleanCategory temp
fi

if [ "$DIRECTIVE" == "onlyclean" ]
then
    cleanCategory temp
    exit
fi

if [ "$DIRECTIVE" == "shinyclean" ]
then
    cleanCategory temp
    cleanCategory result
fi

if [ "$DIRECTIVE" == "onlyshinyclean" ]
then
    cleanCategory temp
    cleanCategory result
    exit
fi

# set desired number of concurrent processes
# prefer an already set value of $NPROC, or use nproc to return it if available

if [ -z "$NPROC" ]
then
    NPROCBIN=`which nproc`
    if [ -x $NPROCBIN ] 
    then
	# linux with modern core-utils
	NPROC=`$NPROCBIN`
    elif [ -x /proc/cpuinfo ]
    then
	# linux with a proc fs
	NPROC=`grep -c processor /proc/cpuinfo` 
    elif [ -x /usr/sbin/sysctl ]
    then
        # Mac OS X
	NPROC=`/usr/sbin/sysctl -n hw.ncpu`
    fi
    
    if [ -z "$NPROC" ] 
    then 
	NPROC=1
    fi
fi
