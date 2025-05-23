FROM debian:12.10

RUN apt-get update && apt-get -y install python3 python3-pip sickle fastqc ncbi-blast+

RUN mkdir -p /usr/src/tnatlas/src
COPY ./pyproject.toml ./MANIFEST.in README.md ./testrun.sh /usr/src/tnatlas/
COPY ./src /usr/src/tnatlas/src/

RUN chmod u+x /usr/src/tnatlas/testrun.sh
RUN cp /usr/src/tnatlas/testrun.sh /usr/bin/testrun


RUN python3 -m pip install /usr/src/tnatlas --break-system-packages

RUN mkdir /data
COPY ./data /data

CMD /bin/bash