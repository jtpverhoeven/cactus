import socket
import time
import sys
import json
import struct
import json
import os
import base64
import configparser

class directorObject:

    def __init__(self):
        self.host = None
        self.port = None

    def connect(self, host, port):
        self.host = host
        self.port = port
        self.cactusSocket = socket.socket()
        self.cactusSocket.settimeout(5)
        self.cactusSocket.connect((self.host,self.port))
        self.cactusSocket.settimeout(120)

    def sendCommand(self, cmd):
        cmdJSON = json.dumps({'cmd' : cmd})
        cmdJsonBytes = bytes(cmdJSON, 'UTF-8')
        length = len(cmdJsonBytes)
        self.cactusSocket.sendall(struct.pack('!I', length))
        self.cactusSocket.sendall(cmdJsonBytes)
        response = self.cactusSocket.recv(1024).decode()
        return response

    def sendBirectionalCommand(self, cmd):
        cmdJSON = json.dumps({'cmd' : cmd})
        cmdJsonBytes = bytes(cmdJSON, 'UTF-8')
        length = len(cmdJsonBytes)
        self.cactusSocket.sendall(struct.pack('!I', length))
        self.cactusSocket.sendall(cmdJsonBytes)

        if cmd != 'bye':
            response = self.cactusSocket.recv(1024).decode()

        return response

    def sendFile(self, sendFile, destFile):
        with open(sendFile, 'r') as f:
            dataToSend = f.read()
            f.close()

        cmdJSON = json.dumps({'cmd' : 'receive_file',  'destination_file': destFile, 'contents': dataToSend})
        cmdJsonBytes = bytes(cmdJSON, 'UTF-8')
        length = len(cmdJsonBytes)
        self.cactusSocket.sendall(struct.pack('!I', length))
        self.cactusSocket.sendall(cmdJsonBytes)
        response = self.cactusSocket.recv(1024).decode()

    def receiveResults(self):
        cmdJSON = json.dumps({'cmd' : 'return_results'})
        cmdJsonBytes = bytes(cmdJSON, 'UTF-8')
        length = len(cmdJsonBytes)

        self.cactusSocket.sendall(struct.pack('!I', length))
        self.cactusSocket.sendall(cmdJsonBytes)

        data = self.recvall(self.cactusSocket, 4)
        length, = struct.unpack('!I', data)
        dataReceived = self.recvall(self.cactusSocket, length)
        dataObject = json.loads(dataReceived.decode())
        return dataObject



    def recvall(self, sock, count):
        buf = b''
        bufText = ''
        while count:
            newbuf = sock.recv(count)
            if not newbuf: return None
            buf += newbuf
            count -= len(newbuf)
        return buf

    def close(self):
        self.cactusSocket.close()
