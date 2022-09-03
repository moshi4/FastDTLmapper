FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y python3 python3-pip python2 libgomp1 qt5-default

# Install FastDTLmapper & Clear dependencies cache
RUN pip install -U pip && \
    pip install fastdtlmapper --no-cache-dir

CMD ["/bin/bash"]
