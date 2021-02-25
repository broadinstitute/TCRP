FROM conda/miniconda2
RUN conda install pytorch=0.4.1 -c pytorch
RUN pip install scipy==0.16
RUN conda install scikit-learn=0.20 -c conda-forge
RUN conda install h5py
RUN conda install pandas
COPY ./code code
COPY ./reproduce reproduce
