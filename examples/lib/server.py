import BaseHTTPServer, SimpleHTTPServer
import ssl

httpd = BaseHTTPServer.HTTPServer(('brassica.ucsd.edu', 4443), SimpleHTTPServer.SimpleHTTPRequestHandler)
httpd.socket = ssl.wrap_socket (httpd.socket, certfile='/Users/sfederow/ome/om/examples/lib/mycert.pem', server_side=True)
httpd.serve_forever()
