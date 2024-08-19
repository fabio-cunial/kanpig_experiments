FROM continuumio/miniconda3:24.1.2-0
MAINTAINER Fabio Cunial
ARG work_dir=/kanpig_experiments
WORKDIR ${work_dir}


# OS
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        build-essential \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        python2 \
        python-dev \
        bsdmainutils \
        default-jre \
        default-jdk \
        unzip \
        python3-pycurl \
        python3-edlib 
        

# GSUTIL
ARG gcloud_version=405.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
ARG htslib_version=1.19.1
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
ARG samtools_version=1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
ARG bcftools_version=1.19
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && bcftools --help

# SNIFFLES
ARG sniffles_version=2.3.3
RUN pip3 install sniffles==${sniffles_version} \
    && sniffles -h

# CUTESV
RUN pip install setuptools==58 scikit-learn \
    && git clone https://github.com/tjiangHIT/cuteSV.git \
    && cd cuteSV \
    && python setup.py install \
    && cuteSV -h

# KANPIG
ARG kanpig_version=0.3.1
RUN wget https://github.com/ACEnglish/kanpig/releases/download/v0.3.1/kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && tar -xzvf kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && mv kanpig-v${kanpig_version}-x86_64-unknown-linux-musl/kanpig ./kanpig \
    && chmod +x ./kanpig \
    && ./kanpig -h

# LRCALLER 1.0
RUN wget https://github.com/DecodeGenetics/LRcaller/releases/download/v1.0/lrcaller-1.0_linux-x86-64.tgz \
    && tar -xzvf lrcaller-1.0_linux-x86-64.tgz \
    && chmod +x lrcaller-1.0/bin/* \
    && mv lrcaller-1.0/bin/* /usr/local/bin \
    && lrcaller --help || echo "LRcaller returns an error code when called"

# TRUVARI
ARG truvari_version=4.2.2
RUN pip3 install truvari==${truvari_version} \
    && truvari --help

# BEDTOOLS
ARG bedtools_version=2.31.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools.static \
    && mv bedtools.static bedtools \
    && chmod +x bedtools \
    && mv bedtools /usr/local/bin \
    && bedtools --help

# MINIGRAPH
ARG minigraph_version=0.21
RUN wget https://github.com/lh3/minigraph/releases/download/v0.21/minigraph-${minigraph_version}_x64-linux.tar.bz2 \
    && tar xjf minigraph-${minigraph_version}_x64-linux.tar.bz2 \
    && rm -f minigraph-${minigraph_version}_x64-linux.tar.bz2 \
    && chmod +x minigraph-${minigraph_version}_x64-linux/minigraph
ENV PATH=${work_dir}/minigraph-${minigraph_version}_x64-linux:${PATH}
RUN minigraph && echo 0 || echo 1

# SVJEDI-GRAPH 1.2.1 with Adam's patches to work around the need for END
ADD svjedigraph_patched ${work_dir}/svjedigraph

# MONITORING SCRIPT
COPY ./vm_local_monitoring_script.sh .
COPY ../manuscript/full_analysis/analysis.py .