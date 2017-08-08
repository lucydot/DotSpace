from flask import Flask
app = Flask(__name__)

@app.route('/')
def index():
    return 'Index Page'

@app.route('/hello')
def hello():
    author = "Me"
    name = "You"
    return render_template('index.html', author=author, name=name)

if __name__ = "__main__":
    app.run()
