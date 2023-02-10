setup:
	@conda env create --file environment.yml || conda env update --file environment.yml

setup_platon_database:
	@mkdir -p data/dbs/platon 
	@wget -O data/dbs/platon/platon.tar.gz https://zenodo.org/record/4066768/files/db.tar.gz && \
	cd data/dbs/platon && \
	tar -xvf platon.tar.gz && \
	rm platon.tar.gz

setup_abricate_database:
	@abricate --setupdb

setup_bacmet2_database:
	@mkdir -p data/dbs/bacmet2
	@wget -O data/dbs/bacmet2/BacMet2_EXP_database.fasta http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta
	@wget -O data/dbs/bacmet2/BacMet2_predicted_database.fasta.gz http://bacmet.biomedicine.gu.se/download/BacMet2_predicted_database.fasta.gz && \
		cd data/dbs/bacmet2/ && \
		gzip -d -f BacMet2_predicted_database.fasta.gz
	@wget -O data/dbs/bacmet2/BacMet2_EXP.753.mapping.txt http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.753.mapping.txt
	@wget -O data/dbs/bacmet2/BacMet2_PRE.155512.mapping.txt http://bacmet.biomedicine.gu.se/download/BacMet2_PRE.155512.mapping.txt
	@makeblastdb -in data/dbs/bacmet2/BacMet2_EXP_database.fasta -dbtype prot -out data/dbs/bacmet2/bacmet2_exp
	@makeblastdb -in data/dbs/bacmet2/BacMet2_predicted_database.fasta -dbtype prot -out data/dbs/bacmet2/bacmet2_pred

setup_databases: setup_platon_database setup_bacmet2_database setup_abricate_database
