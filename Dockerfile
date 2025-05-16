FROM debian:12.10

RUN apt-get update && apt-get -y install python3 python3-pip sickle fastqc ncbi-blast+

RUN mkdir -p /usr/src/tnatlas/src
COPY ./pyproject.toml ./MANIFEST.in README.md /usr/src/tnatlas/
COPY ./src /usr/src/tnatlas/src/

RUN python3 -m pip install /usr/src/tnatlas --break-system-packages

RUN mkdir /data
COPY ./data /data

CMD /bin/bash