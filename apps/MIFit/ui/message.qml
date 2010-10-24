import Qt 4.7

Item {
    id: messageArea
    width: text.width + 10
    height: text.height + 5
    x: Math.ceil((root.width - width)/2.0)
    anchors.leftMargin: 5
    anchors.bottom: root.bottom
    anchors.bottomMargin: 10

    Text {
        id: text
        anchors.centerIn: parent
        color: "white"
        text: message
    }
    Rectangle {
        id: messageBox
        color: "gray"
        radius:  5
        opacity: 0.5
        anchors.fill: parent
        z: -1
    }

    states: [
        State {
            name: "blank"
            when:  !messageVisible
            PropertyChanges {
                target: messageArea
                opacity: 0
            }
        }
    ]

    transitions: [
        Transition {
            from: ""; to: "blank"; reversible: true
            PropertyAnimation { target: messageArea; property: "opacity"; duration: 300 }
        }
    ]
}
