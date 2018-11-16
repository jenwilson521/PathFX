# Installing UMLS

The following procedures are necessary:

## Read the Perl Module INSTALL First

Download `UMLS::Interface` and `UMLS::Similarity` from MetaCPAN

- https://metacpan.org/pod/UMLS::Similarity
- https://metacpan.org/pod/UMLS::Interface

Unpack `UMLS-Similarity.tar.gz` and read the `INSTALL` instructions carefully.

## Register and Download UMLS

It takes a few days for the UMLS support team to authorize your registration and give you access to the data.

### Run the Java-based Selector

Once you've downloaded the 2017AA UMLS version, you'll need to choose which subsets you want to install. Level 0 requires no additional licenses and installs by default. Start with Level 0 by typing:

	sh run_mac.sh

Follow the GUI prompts. Be sure to select **install MySQL loaders**; when everything is configured, go to the **Done** menu and start installing. Note: this step is not mentioned anywhere in the UMLS install instructions, but fortunately it *is* mentioned in the `UMLS::Interface` install instructions.

### Create a Database and Anonymous User

The anonymous user must be configured to have no password. This is a really bad idea in terms of security, but it's what `UMLS::Interface` wants. Lastly, the user can be contained to only access `umls` and `umlsinterfaceindex` databases.

	mysql -u root
	> create database umls;
	> grant all privileges on umls.* to ''@'%';

### Create a Database for UMLS::Interface Testing

Nowhere in the documentation is this mentioned, but it's in the
tests.

	mysql -u root
	> create database umlsinterfaceindex;
	> grant all privileges on umlsinterfaceindex.* to ''@'%';


### Load the Data from the MySQL Loaders

You'll need to repeat the next three steps below this for the `NET` and `META` directories.

#### Enable Loading from a Local Infile

Modify the parameters in `populate_net_mysql_db.sh` to appear as below

	MYSQL_HOME=/usr/local
	db_name=umls

Modify the `mysql` command in `populate_net_mysql_db.sh` to appear as below

	$MYSQL_HOME/bin/mysql -vvv --local-infile $db_name < mysql_net_tables.sql >> mysql_net.log 2>&1

#### Update the SQL Commands

Evidently the MySQL loader is out-of-date with current MySQL definition of SQL. If so, this is MySQL's fault for API changes without backwards compatibility.

Modify `mysql_net_tables.sql` with the following commands

	perl -pi -e 's/\t/ /g' mysql_net_tables.sql

- Tabs are not allowed in MySQL v5.7+, released Jan 2018

#### Run the MySQL Loaders

	sh populate_net_mysql_db.sh
	
As mentioned before, you'll need to repeat the previous three steps for both the `NET` and `META` directories.

## Install the Perl Libraries

### UMLS::Interface

Follow the steps manually for best control, just like the good ol' days :)

	perl Makefile.PL
	make
	make test
	sudo make install
	
Note that the `make test` command may take a few hours.

### UMLS::Similarity

Same as above, you may have to `cpanm` some prerequisites.

	perl Makefile.PL
	make
	make test
	sudo make install
	
Also note that the `make test` command may also several hours.


