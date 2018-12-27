from flask import render_template, url_for
from app import app

import pandas as pd
import json
import re

@app.route('/')
@app.route('/index')
def index():

	data = pd.read_csv("app/static/data/bio_sra.csv", delimiter='\t')
	# data = pd.read_csv(url_for('static', filename='data/bio_sra.csv'), delimiter='\t')
	COLUMNS = ["Bioproject", "species", "env_feature", "env_biome"]
	subdata = data[COLUMNS]
	subdata = subdata[(subdata["species"].notnull()) & 
	                  (subdata["env_feature"].notnull())].reset_index(drop=True)


	subdata["species"] = subdata["species"].apply(lambda species: re.sub('[^A-Za-z0-9]+', '', species))
	print(subdata.shape)
	# out_data = subdata.head().to_json(orient="records")
	out_data = subdata.head(1000).to_json(orient="records")
	# out_data = subdata.to_json(orient="records")
	# out_data = subdata.head(2)
	return render_template('index.html', title='Home',  data=out_data)





