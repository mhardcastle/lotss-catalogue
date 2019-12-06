from panoptes_client import Subject, SubjectSet, Panoptes
from panoptes_client.panoptes import PanoptesAPIException
import getpass
import matplotlib.pyplot as mplplot
import pandas as pd
import argparse


def get_image_height(imagePath):
    return mplplot.imread(imagePath).shape[1]


def make_metadata(ellipses, imagePath):
    metaData = {}
    for ellipseCounter, ellipse in ellipses.iterrows():
        metaData.update(
            {
                f"#feedback_{ellipseCounter}_id": ellipseCounter + 1,
                f"#feedback_{ellipseCounter}_x": ellipse["centre_x"],
                f"#feedback_{ellipseCounter}_y": get_image_height(imagePath)
                - ellipse["centre_y"],
                f"#feedback_{ellipseCounter}_toleranceA": ellipse["a"],
                f"#feedback_{ellipseCounter}_toleranceB": ellipse["b"],
                f"#feedback_{ellipseCounter}_theta": ellipse["theta"],
            }
        )
    return metaData


def make_tutorial_images(imagePaths, ellipseData, projectData):
    # Connect to Panoptes
    Panoptes.connect(
        username=projectData["user_name"], password=projectData["password"]
    )

    newSubjects = []
    for imageId, imagePath in enumerate(imagePaths):
        print(f"Adding {imagePath}...")
        try:
            subjectSet = SubjectSet.find(projectData["subject_set"])
        except PanoptesAPIException as e:
            print(e)
            return
        newSubject = Subject()
        newSubject.add_location(imagePath)
        newSubject.links.project = subjectSet.links.project
        newSubject.metadata.update(
            make_metadata(
                ellipseData.get_group(imageId).reset_index(drop=True), imagePath
            )
        )
        newSubject.save()
        newSubjects.append(newSubject)
    subjectSet.add(newSubjects)


if __name__ == "__main__":
    print("Launching traing image uploader...")

    parser = argparse.ArgumentParser(
        description="Upload tutorial subjects with feedback metadata"
    )

    parser.add_argument(
        "--imagePaths",
        "-i",
        metavar="PATH",
        required=True,
        type=str,
        action="store",
        help="path to text file containing list of subject image paths (one path per line)",
    )

    parser.add_argument(
        "--ellipseData",
        "-e",
        metavar="PATH",
        required=True,
        type=str,
        action="store",
        help="path to CSV path containing ellipse data for all images. Required columns: image_id, centre_x, centre_y, a, b, theta",
    )

    parser.add_argument(
        "--subjectSetId",
        "-s",
        metavar="ID",
        required=True,
        type=int,
        action="store",
        help="Zooniverse subject set ID. You must create the subject set using the Project Builder.",
    )

    parser.add_argument(
        "--user",
        "-u",
        metavar="USERNAME",
        type=str,
        action="store",
        help="Zooniverse username.",
    )

    parser.add_argument(
        "--password",
        "-p",
        metavar="PASSWORD",
        type=str,
        action="store",
        help="Zooniverse password.",
    )

    args = parser.parse_args()

    username = input(f"Enter Zooiverse username:") if args.user is None else args.user
    password = (
        getpass.getpass(f"Enter Zooniverse password for {username}: ")
        if args.password is None
        else args.password
    )

    projectData = {
        "subject_set": args.subjectSetId,
        "user_name": username,
        "password": password,
    }

    imagePaths = pd.read_csv(args.imagePaths, header=None)[0]
    ellipseData = pd.read_csv(args.ellipseData).groupby(by="image_id")

    try:
        print(
            "Adding subjects to {} {}".format(
                SubjectSet.find(projectData["subject_set"]).display_name,
                projectData["subject_set"],
            )
        )
    except PanoptesAPIException as e:
        print(e)
        exit(1)

    make_tutorial_images(imagePaths, ellipseData, projectData)
