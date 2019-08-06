from flask import Flask, render_template

from . import db, run, download

app = Flask(__name__)
@app.route('/')
def index():
    return render_template('index.html')

app.register_blueprint(db.bp)
app.register_blueprint(run.bp)
app.register_blueprint(download.bp)
