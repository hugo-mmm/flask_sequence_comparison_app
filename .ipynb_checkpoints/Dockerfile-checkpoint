# Use an official Python runtime as a parent image
FROM python:3.10

# Set the working directory in the container to /app
WORKDIR /app

# Add the current directory contents into the container at /app
ADD . /app

# Set up a virtual environment
RUN python -m venv venv

# Install any needed packages specified in requirements.txt
RUN /bin/bash -c "source venv/bin/activate && pip install --no-cache-dir -r requirements.txt"

# Make port 80 available to the world outside this container
EXPOSE 80

# Run app.py when the container launches
CMD [ "venv/bin/python", "app.py" ]