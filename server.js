var budo = require('budo')
var babelify = require('babelify')

budo('./src/index.js', {
  host: "localhost",
  port: 8000,             // use this port
  browserify: {
    transform: babelify   // ES6
  }
}).on('connect', function (ev) {
  console.log('Server running on %s', ev.uri)
}).on('update', function (buffer) {
  console.log('bundle - %d bytes', buffer.length)
})
