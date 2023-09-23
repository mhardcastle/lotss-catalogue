#!/usr/bin/python
# Generic file download with retry and check for length

from __future__ import print_function
import requests
import os
from time import sleep,time

def download_file(url,filename):
    if os.path.isfile(filename):
        print('File',filename,'already exists, skipping')
        return
    else:
        print('Downloading',filename)

    downloaded=False
    while not downloaded:
        connected=False
        while not connected:
            try:
                response = requests.get(url, stream=True,timeout=60,verify=False)
                if response.status_code>=502 and response.status_code<=504:
                    print('Code %i -- retrying in 10 seconds!' % response.status_code)
                    sleep(10)
                    continue
                if response.status_code!=200:
                    raise RuntimeError('Code was %i' % response.status_code)
                try:
                    esize=int(response.headers['Content-Length'])
                except KeyError:
                    esize=None
            except (requests.exceptions.ConnectionError,requests.exceptions.Timeout,requests.exceptions.ReadTimeout):
                print('Connection error! sleeping 30 seconds before retry...')
                sleep(30)
            else:
                connected=True
        try:
            if esize is not None:
                print('Downloading %i bytes' % esize)
            else:
                print('Downloading file of unknown length')
            starttime=time()
            with open(filename, 'wb') as fd:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        fd.write(chunk)
            fsize=os.path.getsize(filename)
            if esize!=fsize and esize is not None:
                print('Download incomplete (expected %i, got %i)! Retrying' % (esize, fsize))
            else:
                endtime=time()
                dt=endtime-starttime
                print('Download successful, %i of %s bytes received in %.2f seconds (%.2f MB/s)' % (fsize, str(esize) if esize is not None else 'unknown', dt, fsize/(dt*1024*1024)))
                downloaded=True

        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ChunkedEncodingError):
            print('Connection error! sleeping 30 seconds before retry...')
            sleep(30) # back to the connection


    del response
    return downloaded
