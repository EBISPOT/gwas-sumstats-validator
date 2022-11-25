FROM python:3.7-slim-buster

RUN apt-get update \
    && apt-get install -y --no-install-recommends gcc python-dev \
    && rm -rf /var/lib/apt/lists/*

COPY setup.py .
COPY requirements.txt .
COPY ss_validate ss_validate
COPY tests tests
COPY README.md .


RUN pip install -r requirements.txt \
    && pip install . \
    && apt-get purge -y --auto-remove gcc python-dev