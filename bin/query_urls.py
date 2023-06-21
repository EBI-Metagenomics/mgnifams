#!/usr/bin/env python3

import requests
import sys

urls_file = sys.argv[1]

with open(urls_file) as urls:
    for url in urls:
        url = url.strip()
        response = requests.get(url)
        protein_id = url.split('/')[-1]
        with open(f"{protein_id}_page.html", 'w') as outfile:
            outfile.write(response.text)
