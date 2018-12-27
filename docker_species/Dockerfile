FROM amancevice/pandas:0.23.4 


WORKDIR /home/app
COPY . .
RUN pip install -r requirements.txt

RUN chmod +x boot.sh


EXPOSE 5000
ENTRYPOINT ["./boot.sh"]